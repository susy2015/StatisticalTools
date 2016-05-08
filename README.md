# StatisticalTools
==> To produce the limits using the Higgs combination tool

] Following the https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit#ROOT5_SLC6_release_CMSSW_7_1_X to checkout the Higgs combination tool

] git clone git@github.com:susy2015/StatisticalTools.git

] cd StatisticalTools/CombCardsMaker/python

] run an example of T2tt

  cd batch_T2tt_37Bins/flatNtp_v6/

  ln -s /uscms_data/d3/lhx/publik/for_StatisticalTools/v3_add_data_vs_mc_tagger_unc signalCards

  ./batchCombJobs.py -s signalCards -o prod >& log_prod.lg&

  ./onlyMakeCards.py -s signalCards

==> To draw the limit plots

cd PlotsSMS/workdir/

make -j9

cd T2tt_37Bins/flatNtp_v6/

source batch.txt
