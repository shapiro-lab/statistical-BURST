# statistical-BURST
Public repository for scripts to statistically process BURST images

## Data import

The data immport section currently included in the code assumes two types of input data: intensity images (`is_IQ = false`) or beamformed IQ images (`is_IQ = true`). After setting the data type (`is_IQ`), the script reads all `.mat` files from a specified folder and assembles them into an image sequence (`image_seq`). Folder structure follows:

```
your folder/
  ├── image_block_1.mat
  ├── image_block_2.mat
  ├── image_block_3.mat
  ├── ...
  ├── IQ_block_1.mat (optional)
  ├── IQ_block_2.mat (optional)
  ├── IQ_block_3.mat (optional)
  └── ...
```

Please set 

`basedir = 'enter your folder';`

Each `.mat` file should contain one frame of data. For intensity images, each file is assumed to contain `RData`. For IQ data, each file is assumed to contain `IQ`, whose real and complex parts are stacked. You can also modify this data import section based on your folder structure.
