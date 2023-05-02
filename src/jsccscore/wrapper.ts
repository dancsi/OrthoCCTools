declare function callMain(args: string[]): number;

var Module = {
    "print": (text: string) => { console.log(text) },
    "printErr": (text: string) => { console.log(text) },
    "noInitialRun": true,
    "onRuntimeInitialized": () => {
        console.log("runtime initialized");
        initializeUI();
    }
};

enum Orientation {
    Parallel = 1,
    Antiparallel = 2,
    Both = Parallel | Antiparallel
};

enum ScoreFunctions {
    Potapov = "potapov",
    bCIPA = "bcipa",
    qCIPA = "qCIPA",
    iCIPA_core_vert = "icipa_core_vert",
    iCIPA_nter_core = "icipa_nter_core"
};

function getBaseName(fileName: string) {
    const lastDot = fileName.lastIndexOf(".");
    return fileName.slice(0, lastDot);
}

interface FastscoreArguments {
    inputFile: File,
    maxHeptadDisplacement: number,
    truncate: boolean,
    orientation: Orientation,
    scoreFunction: ScoreFunctions
}

async function fastscoreCommandlineWrapper(
    { 
        inputFile,
        maxHeptadDisplacement = 0,
        truncate = false,
        orientation = Orientation.Parallel,
        scoreFunction = ScoreFunctions.Potapov 
    }: FastscoreArguments
): Promise<File[]> {

    const inputBufView = new Uint8Array(await inputFile.arrayBuffer());
    FS.writeFile(inputFile.name, inputBufView);

    const inputFileName = inputFile.name;
    const baseName = getBaseName(inputFileName);

    const argv = [inputFile.name, `--max-heptad-displacement=${maxHeptadDisplacement}`, `--truncate=${truncate ? 1 : 0}`, `--orientation=${Orientation[orientation]}`, `--score-func=${scoreFunction}`];
    console.log(argv);

    callMain(argv);

    const generatedFiles: File[] = [];

    for (const path of FS.readdir(".")) {
        if (path.startsWith(baseName) && path != inputFileName) {
            generatedFiles.push(new File([FS.readFile(path)], path));
        }
    }

    return generatedFiles;
}