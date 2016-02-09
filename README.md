# StatisticalTools
==> To produce the limits using the Higgs combination tool

] Following the https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit#ROOT5_SLC6_release_CMSSW_7_1_X to checkout the Higgs combination tool
] git clone git@github.com:susy2015/StatisticalTools.git
] cd StatisticalTools/CombCardsMaker/python
] make the several links:

ln -s datacard_v3/DataCard_QCD_v4.txt qcd.txt
ln -s datacard_v3/data.txt .
ln -s datacard_v3/ttz.txt .
ln -s datacard_v3/HadTau_card_RelSys.txt hadtau.txt
ln -s datacard_v3/zinv_datacard.txt zinv.txt
ln -s datacard_v3/LL_card.txt lostle.txt

ln -s signalScan_prod/v3 signalScan
ln -s ../../PlotsSMS/Utils/xSec_T3G.root xSec_T3G.root

] submit the jobs (all the jobs are done in the batchCombJobs.py. We'll get final results to be used to draw the limits.
./batchCombJobs.py -s signalScan -o prod

==> To draw the limit plots
cd PlotsSMS/workdir/
ln -s ../Utils .
ln -s ../python .
rm ori_combined_plots_T2tt
ln -s ../../CombCardsMaker/python/cards_for_plotting ori_combined_plots_T2tt
ln -s Utils/* .
make
source batch.txt
