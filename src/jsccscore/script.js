function encodeText(memory, base, string) {
    for (let i = 0; i < string.length; i++) {
        memory[base + i] = string.charCodeAt(i);
    }
    memory[base + string.length] = 0;
}

function decodeText(memory, base) {
    let cursor = base;
    let result = '';

    while (memory[cursor] !== 0) {
        console.log(memory[cursor])
        result += String.fromCharCode(memory[cursor++]);
    }

    return result;
}

function prepareStrings(memory, buf1, buf2, seq1, seq2) {
    encodeText(memory, buf1, seq1);
    encodeText(memory, buf2, seq2);
}

function setChangeCallback(callback) {
    const elSeq1 = document.getElementById("seq1");
    const elSeq2 = document.getElementById("seq2");
    const elAlign = document.getElementById("align");
    const elTruncate = document.getElementById("truncate");
    let orientation = document.querySelector('input[name="orientation"]:checked').id
    let scoringEngine = document.querySelector('input[name="scoringEngine"]:checked').id;

    const elRes = document.getElementById("res");


    const realCallback = () => {
        const seq1 = elSeq1.value;
        const seq2 = elSeq2.value;
        const align = elAlign.valueAsNumber;
        const truncate = elTruncate.checked ? 1 : 0;

        const res = callback(seq1, seq2, align, truncate, orientation, scoringEngine);
        elRes.textContent = res.toString();
    }

    realCallback();

    elSeq1.onchange = realCallback;
    elSeq2.onchange = realCallback;
    elAlign.onchange = realCallback;
    elTruncate.onchange = realCallback;
    document.getElementById("orientationChooser").addEventListener("change", (e) => {
        orientation = e.target.id;
        realCallback();
    }, false);
    document.getElementById("scoringEngineChooser").addEventListener("change", (e) => {
        scoringEngine = e.target.id;
        realCallback();
    }, false);
}

function main() {
    const orientationLookup = {
        "parallel": 0,
        "antiparallel": 1,
        "both": 2
    };

    const memory = new WebAssembly.Memory({ initial: 2 });
    (async () => {
        const wasmBuf = await (await fetch("libscore.wasm")).arrayBuffer();
        const { instance } = await WebAssembly.instantiate(
            wasmBuf, { "env": { "memory": memory } }
        );
        let maxLen = 1024;
        let buf1 = instance.exports.alloc_string(maxLen);
        let buf2 = instance.exports.alloc_string(maxLen);   

        setChangeCallback((seq1, seq2, align, truncate, orientationStr, engineStr) => {
            const len = Math.max(seq1.length, seq2.length);
            if (len > maxLen) {
                instance.exports.free_string(buf1);
                instance.exports.free_string(buf2);
                maxLen = len;
                buf1 = instance.exports.alloc_string(maxLen);
                buf2 = instance.exports.alloc_string(maxLen);
            }
            const view = new Uint8Array(memory.buffer);
            prepareStrings(view, buf1, buf2, seq1, seq2);

            const orientation = orientationLookup[orientationStr];
            const engine = instance.exports[engineStr];
            return engine(buf1, maxLen, buf2, maxLen, align, truncate, orientation);
        });
    })();
}

main()
