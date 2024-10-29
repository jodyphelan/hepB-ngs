# hepB-ngs

## Install 

```
wget https://raw.githubusercontent.com/jodyphelan/hepB-ngs/main/conda/linux-env.txt
conda create --name hepB --file linux-env.txt
conda activate hepB
pip install --force-reinstall  git+https://github.com/jodyphelan/pathogen-profiler.git
pip install --force-reinstall  git+https://github.com/jodyphelan/hepB-ngs.git
```

## Set up reference database
```
hepB-download-ref.py
```

## Run
```
hepB-ngs.py -f folder_with_fastq_files
```

