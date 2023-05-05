var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
var Module = {
    "print": (text) => { console.log(text); },
    "printErr": (text) => { console.log(text); },
    "noInitialRun": true,
    "onRuntimeInitialized": () => {
        console.log("runtime initialized");
        initializeUI();
    }
};
var Orientation;
(function (Orientation) {
    Orientation[Orientation["Parallel"] = 1] = "Parallel";
    Orientation[Orientation["Antiparallel"] = 2] = "Antiparallel";
    Orientation[Orientation["Both"] = 3] = "Both";
})(Orientation || (Orientation = {}));
;
var ScoreFunctions;
(function (ScoreFunctions) {
    ScoreFunctions["Potapov"] = "potapov";
    ScoreFunctions["bCIPA"] = "bcipa";
    ScoreFunctions["qCIPA"] = "qCIPA";
    ScoreFunctions["iCIPA_core_vert"] = "icipa_core_vert";
    ScoreFunctions["iCIPA_nter_core"] = "icipa_nter_core";
})(ScoreFunctions || (ScoreFunctions = {}));
;
var OutputFormat;
(function (OutputFormat) {
    OutputFormat["Binary"] = "bin";
    OutputFormat["CSV"] = "csv";
})(OutputFormat || (OutputFormat = {}));
;
function getBaseName(fileName) {
    const lastDot = fileName.lastIndexOf(".");
    return fileName.slice(0, lastDot);
}
function fastscoreCommandlineWrapper({ inputFile, maxHeptadDisplacement = 0, truncate = false, orientation = Orientation.Parallel, scoreFunction = ScoreFunctions.Potapov, outputFormat = OutputFormat.Binary }) {
    return __awaiter(this, void 0, void 0, function* () {
        const inputBufView = new Uint8Array(yield inputFile.arrayBuffer());
        FS.writeFile(inputFile.name, inputBufView);
        const inputFileName = inputFile.name;
        const baseName = getBaseName(inputFileName);
        const argv = [inputFile.name, `--max-heptad-displacement=${maxHeptadDisplacement}`, `--truncate=${truncate ? 1 : 0}`, `--orientation=${Orientation[orientation]}`, `--score-func=${scoreFunction}`, `--output-format=${outputFormat}`];
        showCommandLine(argv);
        console.log(argv);
        callMain(argv);
        const generatedFiles = [];
        for (const path of FS.readdir(".")) {
            if (path.startsWith(baseName) && path != inputFileName) {
                generatedFiles.push(new File([FS.readFile(path)], path));
            }
        }
        return generatedFiles;
    });
}
