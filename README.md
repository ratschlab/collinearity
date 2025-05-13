
## Build from source

```bash
git clone --recursive https://github.com/ratschlab/collinearity.git
cd collinearity

## install sdsl
#pushd external/sdsl-lite
#./install.sh $PWD
#popd

# for python bindings
pip install .

# for cpp
mkdir build && cd build
cmake ..
make -j4
```

