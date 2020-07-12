# schedule



# files
* dataAnalyzeCriticalAll.mlx preprocess for data_rat010_0615_event, which is used for conference paper. All other data should be preprocessed 
  in the same way.

* *folder* datFiles
  * /dataFiles/data: files named as ratxxx_mmdd_xM1_ymPFC.mat contain spike trains, mapping matrix, segment matrix, etc
  * /dataFiles/results: contain foldes named as ratxxx_mmdd, which contained trained models and tested statistics

* *folder* preprocess
  * each dataset correponds to a file with similar name. The file transfers spike time into spike times and preprocess the spike trians,
    select neruons, ect
* *folder* ratadata
  * raw data files