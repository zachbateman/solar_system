import os
from pycallgraph import PyCallGraph
from pycallgraph import Config
from pycallgraph import GlobbingFilter
from pycallgraph.output import GraphvizOutput

from solar_system import *


config = Config()
config.trace_filter = GlobbingFilter(exclude=['pycallgraph.*',
                                                           '*.secret_function',
                                                           '_*',
                                                           'SourceFileLoader*',
                                                           'ModuleSpec*',
                                                           'ExtensionFileLoader*',
                                                           'FileFinder*',
                                                           'find_spec*',
                                                           '<listcomp>',
                                                           'spec_*',
                                                           'module_*',
                                                           'cache_*',
                                                           '<genexpr>',
                                                           'cb',
                                                           '*FileFinder',
                                                           '<setcomp>',
                                                           ])
graphviz = GraphvizOutput(output_file='PyCallGraph.png')
with PyCallGraph(output=graphviz, config=config):
    x = SolarSim(TOTAL_TIME=3000, INITIAL_SPACEBODIES=150, ANIMATION_INTERVAL=50, FRAME_SAMPLE_RATE=15)
