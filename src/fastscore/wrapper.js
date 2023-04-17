Module["noInitialRun"] = true;
Module["onRuntimeInitialized"] = _ => {
    FS.writeFile("input.fasta", String.raw`>P2
AEIAALEAEIAALEAEIAALEAEIAALK
>P99
AEIAALEAEIAALKAKIAALEAENAALE`);
    callMain(["input.fasta"]);
    console.log(FS.readFile("input.bin"));
};