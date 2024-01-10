import subprocess 

#com1='export GRAFPATH=$CONDA_PREFIX/share/grafpop'
#com0='echo $GRAFPATH'
com2='grafpop /home/liaoks/GwasQcPipeline/with_modifications_run/sample_level/call_rate_2/samples /home/liaoks/GwasQcPipeline/with_modifications_run/sample_level/call_rate_2/testing.txt'

#subprocess.run(com1, shell=True)
#subprocess.run(com0, shell=True)
subprocess.run(com2, shell=True)

# needed at first line of Plot... #!/home/liaoks/miniconda3/envs/grafpop/bin/perl
# NEED TO FIX BUILD.SH TO INCLUDE THIS 

com3='~/software/GrafPop1.0/PlotGrafPopResults.pl ~/GwasQcPipeline/with_modifications_run/sample_level/call_rate_2/testing.txt ~/GwasQcPipeline/with_modifications_run/sample_level/call_rate_2/testing.png'
subprocess.run(com3, shell=True)

