#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_MetaProfile
----------------------------------

Tests for `MetaProfile` module.
"""

import os
import unittest
import time
import sys
import pylab as pl
sys.path.append("..")
from MetaProfile import utils as metafun
from MetaProfile import Signal, Window, Profile, MetaProfiler, get_profiler, get_signals
os.system("taskset -p 0xff %d" % os.getpid())


class TestMetaprofile(unittest.TestCase):

    def setUp(self):
        pass

    def test_something(self):
        pass

    def tearDown(self):
        pass

signal_list = ({'name': "Test-1",
                'filepath': "test1.bed"},
               {'name': "Test-2",
                'filepath': "test2.bed"},
               {'name': "Test-3",
                'filepath': "test3.bed"},
                )
# signal_list = ({'name': "Conservation", 'filepath': "./hg19.100way.phyloP100way.bw",
                # 'filetype': "bigwig"},)
windows_list = ({'name': "PolyA-signal",
                 'filepath': "test_windows_small.bed",
                 'pseudocount': 0},)

print "########################################################################################################"
start_time = time.time()
start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
sys.stderr.write("############## Started script on %s ##############\n" % start_date)

pc = get_profiler(signal_list, windows_list)
fig, ax = pl.subplots()
for profile_name, profile in pc.profiles.iteritems():
    print profile_name, profile.get_aggregated_profile("raw", metric='mean').sum()
    profile.plot_line("normalized_to_gene", ax=ax, aggregate=True, label=profile_name, metric='sum')
    # ax.plot(profile.get_aggregated_profile("raw", metric='mean')/profile.get_aggregated_profile("raw", metric='mean').sum(),
    #        label=profile_name)
pl.legend()
pl.savefig("plot.pdf")
# pc.profiles['Conservation on PolyA-signal'].plot_heatmap("normalized_to_gene",
        # sort_by=pc.profiles['Conservation on PolyA-signal'].profile_normalized_to_gene.mean(axis=1),
                                                          # cmap='Reds')
# pl.savefig("heatmap.pdf")

sys.stderr.write("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))
print "########################################################################################################"
# start_time = time.time()
# start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
# sys.stderr.write("############## Started script on %s ##############\n" % start_date)

# pc = get_profiler(signal_list, windows_list, parallel=False)

# sys.stderr.write("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))
# print "########################################################################################################"

# print "########################################################################################################"
# start_time = time.time()
# start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
# sys.stderr.write("############## Started script on %s ##############\n" % start_date)
# signal_list = ({'name': "CPSF-160-60min",
#                 'filepath': "/import/bc2/home/zavolan/gumiennr/Ule/Ule/CLIPzSamples/AnalysisForPaper/time_60/mapped_sequences.all.bed"},
#                {'name': "CPSF-160-120min",
#                 'filepath': "/import/bc2/home/zavolan/gumiennr/Ule/Ule/CLIPzSamples/AnalysisForPaper/time_120/mapped_sequences.all.bed"},
#                {'name': "CPSF-160-30min",
#                 'filepath': "/import/bc2/home/zavolan/gumiennr/Ule/Ule/CLIPzSamples/AnalysisForPaper/time_30/mapped_sequences.all.bed"},
#                {'name': "CPSF-160-240min",
#                 'filepath': "/import/bc2/home/zavolan/gumiennr/Ule/Ule/CLIPzSamples/AnalysisForPaper/time_240/mapped_sequences.all.bed"})

# windows_list = ({'name': "PolyA-signal",
#                  'filepath': "/import/bc2/home/zavolan/gumiennr/Ule/Ule/AdditionalData/PolyAsignal/PolyAsignal_from_asymetric.bed",
#                  'pseudocount': 0},)
# pc = get_profiler(signal_list, windows_list)
# sys.stderr.write("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))
