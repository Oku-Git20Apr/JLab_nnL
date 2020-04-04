#! /bin/sh

days="1119"

H1="./bin/ana_Lambda -f ../run_list/nnlambda/Lambda_small_optH1.list -p param/f1_tuned_Lambda_twc.param -H1 -m ../matrix/matrix.list -r ../rootfiles/mmass/ana_Lambda/Lambda_small_optH1_$days.root -w ../pdf/mmass/ana_Lambda/Lambda_small_optH1_$days.pdf"

H2="./bin/ana_Lambda -f ../run_list/nnlambda/Lambda_small_optH2.list -p param/f1_Lambda_phase2_tuned.param -H1 -m ../matrix/matrix.list -r ../rootfiles/mmass/ana_Lambda/Lambda_small_optH2_$days.root -w ../pdf/mmass/ana_Lambda/Lambda_small_optH2_$days.pdf"


hadd_H="hadd ../rootfiles/mmass/ana_Lambda/Lambda_small_optH_$days.root ../rootfiles/mmass/ana_Lambda/Lambda_small_optH1_$days.root ../rootfiles/mmass/ana_Lambda/Lambda_small_optH2_$days.root"

H3="./bin/ana_Lambda -f ../run_list/nnlambda/Lambda_small_optT.list -p param/f1_Lambda_phase2_tuned.param -H2 -m ../matrix/matrix.list -r ../rootfiles/mmass/ana_Lambda/Lambda_small_optT_$days.root -w ../pdf/mmass/ana_Lambda/Lambda_small_optT_$days.pdf"

T1="./bin/ana_Lambda -f ../run_list/nnlambda/nnL_small_opt1.list -p param/f1_tuned_Lambda_twc.param -T -m ../matrix/matrix.list -r ../rootfiles/mmass/ana_Lambda/nnL_small_opt1_$days.root -w ../pdf/mmass/ana_Lambda/nnL_small_opt1_$days.pdf"

T2="./bin/ana_Lambda -f ../run_list/nnlambda/nnL_small_opt2.list -p param/f1_Lambda_phase2_tuned.param -T -m ../matrix/matrix.list -r ../rootfiles/mmass/ana_Lambda/nnL_small_opt2_$days.root -w ../pdf/mmass/ana_Lambda/nnL_small_opt2_$days.pdf"

T3="./bin/ana_Lambda -f ../run_list/nnlambda/nnL_small_opt3.list -p param/f1_Lambda_phase2_tuned.param -T -m ../matrix/matrix.list -r ../rootfiles/mmass/ana_Lambda/nnL_small_opt3_$days.root -w ../pdf/mmass/ana_Lambda/nnL_small_opt3_$days.pdf"

T4="./bin/ana_Lambda -f ../run_list/nnlambda/nnL_small_opt4.list -p param/f1_Lambda_phase2_tuned.param -T -m ../matrix/matrix.list -r ../rootfiles/mmass/ana_Lambda/nnL_small_opt4_$days.root -w ../pdf/mmass/ana_Lambda/nnL_small_opt4_$days.pdf"

hadd_T="hadd ../rootfiles/mmass/ana_Lambda/nnL_small_opt_$days.root ../rootfiles/mmass/ana_Lambda/nnL_small_opt1_$days.root ../rootfiles/mmass/ana_Lambda/nnL_small_opt2_$days.root 
 ../rootfiles/mmass/ana_Lambda/nnL_small_opt3_$days.root ../rootfiles/mmass/ana_Lambda/nnL_small_opt4_$days.root"


#T2="./bin/ana_Lambda -f ../run_list/nnlambda/nnL_small_234.list -p param/f1_Lambda_phase2_tuned.param -T -m ../matrix/matrix.list -r ../rootfiles/mmass/ana_Lambda/nnL_small_2_$days.root -w ../pdf/mmass/ana_Lambda/nnL_small_2_$days.pdf"
#hadd_T="hadd ../rootfiles/mmass/ana_Lambda/nnL_small_$days.root ../rootfiles/mmass/ana_Lambda/nnL_small_1_$days.root ../rootfiles/mmass/ana_Lambda/nnL_small_2_$days.root"


eval $H1
eval $H2
eval $H3
eval $T1
eval $T2
eval $T3
eval $T4
eval $hadd_H
eval $hadd_T
