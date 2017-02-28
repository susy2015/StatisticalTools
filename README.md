# StatisticalTools
==> To produce the limits using the Higgs combination tool

] Following the https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit#ROOT6_SLC6_release_CMSSW_7_4_X to checkout the Higgs combination tool

] git clone -b Moriond2017 git@github.com:susy2015/StatisticalTools.git

] cd StatisticalTools/CombCardsMaker/python

] run an example of T2tt

  cd batch_T2tt_Moriond_2017/

  ln -s /uscms_data/d3/lhx/tasks/recipes/prepare2017/ana_Summer16_MC_23Sep2016_Data_Feb07_2017_CMSSW_8_0_25/src/SusyAnaTools/Tools/condor/signalScan/signalCards/T2tt/manu/ signalCards
  
  ln -s ../datacard_v13_Moriond_2017/*.txt .
  ln -s ../new_onlyMakeCards.py .
  ln -s ../new_batchCombJobs.py .
  ln -s ../new_combCardPerChannel.py .
  ln -s ../xSec.root .
  
  python new_batchCombJobs.py -s signalCards -o prod >& log_prod.lg&

  python new_onlyMakeCards.py -s signalCards -m T2tt

==> To draw the limit plots 

cd PlotsSMS/workdir/

make -j9

cd batch_T2tt_Moriond_2017

ln -s $CMSSW_BASE/src/StatisticalTools/CombCardsMaker/python/batch_T2tt_Moriond_2017/cards_for_plotting ori_combined_plots_T2tt

source batch.txt
