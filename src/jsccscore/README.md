To compile the code to WebAssembly, you need the wasi sdk (https://github.com/WebAssembly/wasi-sdk).

```
    export CXXFLAGS="-std=c++17 -I../scoring/ -I.. -lc++abi -lc++ -fno-exceptions -fno-threadsafe-statics -Wl,--no-entry -nostartfiles -Wl,--import-memory -Wl,--allow-undefined -O3 -flto -Wl,--lto-O3"
    export CXX=/opt/wasi-sdk/bin/clang
    $CXX $CXXFLAGS jsccscore.cpp ../scoring/ScoringEnginePotapov.cpp -o libscore.wasm
```