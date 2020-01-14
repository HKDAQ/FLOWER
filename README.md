# energetic-bonsai
Energy reconstruction script that works with (hk-)BONSAI

## Compiling
Make sure the following are sourced
* WCSim (`$WCSIMDIR` set)

Then
```bash
export EBONSAIDIR=/path/to/energetic-bonsai
make
```

## Running
Make sure the following are sourced
* WCSim (`$WCSIMDIR` set)
* BONSAI (`$BONSAIDIR` set)
* energetic_bonsai (`$EBONSAIDIR` set)

Then
```bash
export EBONSAIDIR=/path/to/energetic-bonsai
export PATH=$EBONSAIDIR:$PATH
rootebonsai -b -q energetic_bonsai.C+'("/path/to/wcsim/file.root",1,"SuperK")'
```
