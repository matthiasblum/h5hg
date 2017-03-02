#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import os
import sqlite3
import sys
import tempfile
from subprocess import Popen, PIPE

import h5py
import numpy as np
from scipy.sparse import coo_matrix


class Array:
    def __init__(self, filename):
        self.filename = filename
        self.localqc_res = 500
        self.localqc_data = None
        self.wig_res = 0
        self.wig_data = None
        self.wigu_data = None

    def _get(self, chrom, start, end, localqc=True, wig=True, wigu=True):
        self.localqc_res = 500
        self.localqc_data = None
        self.wig_res = 0
        self.wig_data = None
        self.wigu_data = None

        with h5py.File(self.filename, 'r') as fh:
            grp = fh.get(chrom)

            if grp is None:
                sys.stderr.write('{} not found in {}\n'.format(chrom, self.filename))
                return False

            chrom_size = int(grp.attrs['size'])

            if start < 0:
                start = 0
            if end >= chrom_size:
                end = chrom_size - 1

            if localqc:
                dset = grp.get('localqcs')

                # Select data within interval
                self.localqc_data = dset[start // self.localqc_res:(end + self.localqc_res - 1) // self.localqc_res]

            if wig:
                self.wig_res = int(grp.attrs['span'])
                dset = grp.get('wigs')

                # Select data within interval
                data = dset[start // self.wig_res:(end + self.wig_res - 1) // self.wig_res]
                self.wig_data = data[:, 0]
                if wigu:
                    self.wigu_data = data[:, 1]

        return True

    def get(self, chrom, start, end, localqc=True, wig=True, wigu=True, nreps=0, shrink=0):
        if not self._get(chrom, start, end, localqc, wig, wigu):
            return None

        region_size = end - start

        if localqc:
            if nreps:
                # Use 5-replicates
                flag = {0: 0, 1: 1, 2: 2, 3: 4, 4: 8, 5: 16}.get(nreps, 16)

                # Set dispersion to -1 when localqc does not pass the number of replicates
                self.localqc_data['dispersion'][np.bitwise_and(self.localqc_data['flag'], flag) == 0] = -1
            else:
                # Use 1-replicate

                # Set dispersion to -1 when bin is empty (does not contain any read)
                self.localqc_data['dispersion'][self.localqc_data['intensity'] == 0] = -1

            # Keep only the dispersion
            self.localqc_data = self.localqc_data['dispersion']

            if shrink:
                while self.localqc_res * shrink < region_size:
                    if self.localqc_data.size % 2:
                        self.localqc_data = self.localqc_data[:-1]

                    self.localqc_res *= 2
                    self.localqc_data = self.localqc_data.reshape(self.localqc_data.size // 2, 2).max(axis=1)

                self.localqc_data = self.localqc_data.tolist()

        if wig:
            if shrink:
                while self.wig_res * shrink < region_size:
                    if self.wig_data.size % 2:
                        self.wig_data = self.wig_data[:-1]

                        if wigu:
                            self.wigu_data = self.wigu_data[:-1]

                    self.wig_res *= 2
                    self.wig_data = self.wig_data.reshape(self.wig_data.size // 2, 2).max(axis=1)

                    if wigu:
                        self.wigu_data = self.wigu_data.reshape(self.wigu_data.size // 2, 2).max(axis=1)

            self.wig_data = self.wig_data.tolist()

            if wigu:
                self.wigu_data = self.wigu_data.tolist()

        return {
            'localqc': self.localqc_data,
            'localqcres': self.localqc_res,
            'wig': self.wig_data,
            'wigu': self.wigu_data,
            'wigres': self.wig_res
        }


class Matrix:
    def __init__(self, filename, assembly=None, db=None):
        # HDF5 file
        self.filename = filename

        # Assembly and SQLite database for genes
        self.assembly = assembly
        self.db = db

        # Last genomic region retrieved
        self.chrom = None
        self.start = None
        self.end = None

        # Genes from the last genomic region
        self.genes = []
        self.labels = [[]]

        # Data
        self.rows = None
        self.cols = None
        self.values = None

        # Resolution (bp)
        self.resolution = None

        # Size of the output matrix
        self.matsize = None

        # 90th percentile of retrieved values
        self.p90 = None

        # Boolean indicating if we return contact counts or dispersion
        self.usecounts = True

        # Matrix index of the first bin
        self.x = None

    def get(self, chrom, start=0, end=0, **kwargs):
        sampling = kwargs.get('sampling', 100)
        mincount = kwargs.get('count', 0)
        maxdisp = kwargs.get('disp', 100)
        gene = kwargs.get('gene')

        # If matsize stays None, it means we could not retrieve data
        self.matsize = None

        # Read data from HDF5 file
        with h5py.File(self.filename, 'r') as fh:
            dset = fh.get(chrom)

            if dset is None:
                sys.stderr.write('{} not found in {}\n'.format(chrom, self.filename))
                return self

            chromsize = dset.attrs['size']
            matsize = dset.attrs['matsize']

            # Compute resolution and round it (can be 25002 for 25kb)
            self.resolution = (chromsize - 1) // (matsize - 1)
            self.resolution = self.resolution // 1000 * 1000

            # Update genomic region coordinates
            self.chrom = chrom
            self.start = start if 0 <= start < chromsize else 0
            self.end = end if self.start < end < chromsize else chromsize

            x = self.start // self.resolution
            y = self.end // self.resolution + 1

            self.matsize = y - x

            _dset = dset[
                (dset['row'] >= x) &
                (dset['row'] < y) &
                (dset['col'] >= x) &
                (dset['col'] < y) &
                (dset['data'] >= mincount)
            ]

        try:
            self.p90 = np.percentile(_dset['data'], 90)
        except IndexError:
            self.p90 = 1

        # Filter data on count/dispersion
        if sampling == 100:
            self.usecounts = True

            if maxdisp == 100:
                data = _dset[:]
            else:
                data = _dset[
                    (_dset['disp90'] < maxdisp) &
                    (_dset['disp70'] < maxdisp) &
                    (_dset['disp50'] < maxdisp)
                ]

            self.values = data['data']
        elif sampling == 90:
            self.usecounts = False

            data = _dset[
                (_dset['disp90'] < maxdisp)
            ]

            self.values = data['disp90']
        elif sampling == 70:
            self.usecounts = False

            data = _dset[
                (_dset['disp90'] < maxdisp) &
                (_dset['disp70'] < maxdisp)
                ]

            self.values = data['disp70']
        else:
            self.usecounts = False

            data = _dset[
                (_dset['disp90'] < maxdisp) &
                (_dset['disp70'] < maxdisp) &
                (_dset['disp50'] < maxdisp)
                ]

            self.values = data['disp50']

        self.rows = data['row']
        self.cols = data['col']
        self.x = x

        if self.assembly and self.db:
            self._get_genes(gene, minsize=kwargs.get('minsize', 1), maxlevels=kwargs.get('maxlevels', 5))

        return self

    def _get_genes(self, gene=None, minsize=1, maxlevels=5):
        _genes = []
        self.genes = []
        self.labels = [[]]

        con = sqlite3.connect(self.db)
        cur = con.cursor()
        cur.execute('SELECT symbol, start, end, fw_strand '
                    'FROM gene '
                    'WHERE chrom=:chrom '
                    'AND assembly=:assembly '
                    'AND flag=1 '
                    'AND ((start BETWEEN :min AND :max) '
                    'OR (end BETWEEN :min AND :max) '
                    'OR (:min BETWEEN start AND end)) ', {
                        'chrom': self.chrom,
                        'assembly': self.assembly,
                        'min': self.start,
                        'max': self.end
                    })

        for gene_symbol, tx_start_pos, tx_end_pos, fw_strand in cur:

            if tx_start_pos < self.start:
                tx_start_pos = self.start

            if tx_end_pos > self.end:
                tx_end_pos = self.end

            label_text = gene_symbol + ' >' if fw_strand else '< ' + gene_symbol

            if gene == gene_symbol:
                self.genes.append({
                    'start': tx_start_pos,
                    'end': tx_end_pos,
                    'highlight': True
                })
                self.labels[0].append([tx_start_pos, label_text, True])
            elif tx_end_pos - tx_start_pos >= minsize:
                _genes.append([tx_start_pos, tx_end_pos, label_text])

        cur.close()
        con.close()

        """
        After some tests, I found we can print 117 characters per line
        Characters used: 0-9, a-z, whitespace, "<" symbol
        This will change if the font size or image size are modified
        """
        char_size = float(self.end - self.start) / 117

        # Sort genes by start position
        for tx_start_pos, tx_end_pos, label_text in sorted(_genes, key=lambda x: x[0]):
            self.genes.append({
                'start': tx_start_pos,
                'end': tx_end_pos,
                'highlight': False
            })

            label_start = tx_start_pos
            label_end = tx_start_pos + int(len(label_text) * char_size)

            if label_end >= self.end:
                label_start = tx_end_pos - int(len(label_text) * char_size)

            for i, level in enumerate(self.labels):
                collision = False

                for _label_start, _label_text, _highligh in level:
                    _label_end = _label_start + int(len(_label_text) * char_size)
                    if label_start <= _label_start < label_end \
                            or label_start < _label_end <= label_end \
                            or _label_start < label_end <= _label_end:
                        collision = True
                        break

                if not collision:
                    # Add the label on the current level
                    self.labels[i].append([label_start, label_text, False])
                    break
                elif i + 1 < len(self.labels):
                    # Test on the next level
                    continue
                elif i + 1 < maxlevels:
                    # Add a level
                    self.labels.append([[label_start, label_text, False]])
                    break

    @staticmethod
    def _mkstemp():
        fd, filename = tempfile.mkstemp()
        os.close(fd)
        return filename

    @staticmethod
    def _calc_tics_step(length):
        first_digit = int(str(length)[0])
        if first_digit == 1:
            step = pow(10, math.floor(math.log10(length + 1)) - 1) * 2
        elif first_digit < 5:
            step = pow(10, math.floor(math.log10(length + 1)) - 1) * 5
        else:
            step = pow(10, math.floor(math.log10(length + 1)) - 1) * 10

        return int(step)

    def plot(self, **kwargs):
        astriu = kwargs.get('astriu', False)
        vmax = kwargs.get('vmax')
        width = kwargs.get('width', 1400.)
        colors = kwargs.get('colors',
                            ['#3288BD', '#66C2A5', '#ABDDA4', '#E6F598', '#FEE08B', '#FDAE61', '#F46D43', '#D53E4F'])
        gnuplot = kwargs.get('gnuplot', 'gnuplot')

        if self.matsize is None:
            return None
        elif vmax is None:
            vmax = self.p90

        # Create dense matrix
        if self.usecounts:
            mat_data = coo_matrix(
                (self.values, (self.rows - self.x, self.cols - self.x)),
                shape=(self.matsize, self.matsize),
                dtype=np.int32
            ).todense()
        else:
            mat_data = np.empty(shape=(self.matsize, self.matsize), dtype=np.float32)
            mat_data[np.triu_indices(self.matsize)] = 100  # default value is the highest possible dispersion
            mat_data[self.rows - self.x, self.cols - self.x] = self.values

        # Create a matrix containing the row/col indices
        mat = np.zeros(shape=(self.matsize + 1, self.matsize + 1), dtype=np.float32)
        mat[0, 0] = self.matsize
        mat[0, 1:] = np.arange(self.matsize)
        mat[1:, 0] = np.arange(self.matsize)
        mat[1:, 1:] = np.triu(mat_data) + np.triu(mat_data, 1).T

        # Write the matrix to a temporary binary file
        bin_file = self._mkstemp()
        with open(bin_file, 'wb') as fh:
            fh.write(mat.tostring())

        # Define font/margin settings
        height = width / 2 if astriu else width
        ft_size_default = width * 10 / 1024
        ft_size_cbtics = width * 15 / 1024
        ft_size_labels = width * 8 / 1024
        if astriu:
            ft_size_labels *= 0.75
        lvl_margin = 0.15 if astriu else 0.125

        # Define x-axis tics considering the number of cells to plot on the x-axis
        length = self.end - self.start
        step = self._calc_tics_step(length)
        if length <= self.resolution:
            step *= 2

        tics = []
        for x in range(self.start, self.end + self.resolution + 1, step):
            tics.append((x, float(x - self.start) / self.resolution))

        # Write Gnuplot script
        script = self._mkstemp()
        output = self._mkstemp()
        with open(script, 'w') as fh:
            fh.write("set terminal pngcairo size {},{} enhanced font 'sans,{}'\n".format(width, height, ft_size_default))

            # Output file
            fh.write("set output '{}'\n".format(output))

            # Multiplot
            fh.write('set multiplot\n')

            # Margins
            fh.write('set tmargin 1\n')
            fh.write('set rmargin at screen {}\n'.format(0.85 if astriu else 0.8725))
            fh.write('set bmargin 0\n')
            fh.write('set lmargin at screen 0.05\n')

            # Disable legend
            fh.write('unset key\n')

            # Origin/size of the first plot
            fh.write('set size ratio {}\n'.format(0.5 if astriu else 1))
            fh.write('set origin 0, {}\n'.format(0.1 if astriu else 0.07))

            # Axes/Borders
            fh.write("set style line 11 lc rgb '#808080' lt 1\n")
            fh.write('set border front ls 11\n')
            fh.write('set tics nomirror out scale 0.75\n')

            # No tics for the top plot
            fh.write('unset xtics\n')
            fh.write('unset ytics\n')

            # X/Y-axis limits
            if astriu:
                fh.write('set xrange[-0.5:{}]\n'.format((self.matsize - 1) * 2))
                fh.write("set yrange [-0.5:{}]\n".format(self.matsize - 1))
            else:
                fh.write('set xrange[-0.5:{}]\n'.format(self.matsize - 1))
                fh.write("set yrange [-0.5:{}] reverse\n".format(self.matsize - 1))

            # Color box
            if vmax:
                fh.write("set cbrange [0:{}]\n".format(vmax))

            fh.write('set cbtics font "sans,{}"\n'.format(ft_size_cbtics))

            cblabel = 'contact counts' if self.usecounts else 'dispersion (%)'
            fh.write('set cblabel "{}" offset {},0 font "sans,{}" tc rgb "#808080"\n'.format(
                cblabel,
                2 if astriu else 1,
                ft_size_cbtics
            ))

            # Color palette
            if self.usecounts:
                fh.write("set palette defined ({})\n".format(','.join("{} '{}'".format(i, color) for i, color in enumerate(colors))))
            else:
                fh.write("set palette negative defined ({})\n".format(','.join("{} '{}'".format(i, color) for i, color in enumerate(colors))))

            # Style
            fh.write('set style fill transparent solid 0.8 noborder\n')

            # Plot
            if astriu:
                fh.write("plot '{}' u ($1+$2):($1-$2):3 binary matrix with image\n".format(bin_file))
            else:
                fh.write("plot '{}' binary matrix with image\n".format(bin_file))

            # Margin for 2nd plot
            fh.write('set tmargin 0\n')
            fh.write('set bmargin 0\n')

            # Origin/size
            fh.write('set origin 0, {}\n'.format(0.085 if astriu else 0.05))
            fh.write('set size noratio 1, 0.1\n')

            if astriu:
                fh.write('set object rect from screen 0, screen 0.1 to screen 0.9, screen 0.185 back fs solid noborder\n')
                fh.write('set object rect from screen 0, screen 0 to screen 0.05, screen 1 back fs solid noborder\n')
                fh.write('set object rect from screen 0, screen 0.985 to screen 0.9, screen 1.5 back fs solid noborder\n')

            # Tics
            fh.write('set tics nomirror out scale 0.75\n')
            fh.write('unset ytics\n')
            fh.write('set decimalsign locale "en_US.UTF-8"\n')
            fh.write('set format x "%\'.0f"\n')
            fh.write('set format y "%\'.0f"\n')

            # Limits
            fh.write('set xrange [{}:{}]\n'.format(self.start, self.end))
            fh.write('set yrange[0:1] noreverse\n')

            for gene in self.genes:
                if gene['highlight']:
                    fh.write('set object rect from {},0.8 to {},0.95 fc rgb "#A54F3F" fs solid noborder\n'.format(gene['start'], gene['end']))
                else:
                    fh.write('set object rect from {},0.8 to {},0.95 fc rgb "#46817C" fs solid noborder\n'.format(gene['start'], gene['end']))

            for i, level in enumerate(self.labels):
                y = 0.7 - i * lvl_margin

                for x, label_text, highlight in level:
                    if highlight:
                        fh.write('set label "{}" at {},{} font "sans,{}" tc rgb "#A54F3F"\n'.format(label_text, x, y, ft_size_labels))
                    else:
                        fh.write('set label "{}" at {},{} font "sans,{}"\n'.format(label_text, x, y, ft_size_labels))

            fh.write('plot -1\n')
            fh.write('unset multiplot\n')

        # Run gnuplot with script
        pop = Popen([gnuplot, script], stdout=PIPE, stderr=PIPE)
        out, err = pop.communicate()

        # Delete temporary files
        os.unlink(script)
        os.unlink(bin_file)

        return output

    def todict(self, subtract=False):
        if self.matsize is None:
            return None
        else:
            return {
                'rows': (self.rows - self.x).tolist() if subtract else self.rows.tolist(),
                'cols': (self.cols - self.x).tolist() if subtract else self.cols.tolist(),
                'values': self.values.tolist(),
                'matsize': int(self.matsize),
                'p90': float(self.p90)
            }