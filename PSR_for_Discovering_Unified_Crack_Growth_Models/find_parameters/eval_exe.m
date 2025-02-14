function [eval_out] = eval_exe(evalstr,xtrain_in,Const_pair_now_in)
%EVAL_EXE 
   xtrain=xtrain_in;
   Const_pair_now=Const_pair_now_in;
   eval_out=eval([evalstr ';']);
end

