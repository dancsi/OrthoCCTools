var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
function offerFileForDownload(file) {
    const url = window.URL.createObjectURL(file);
    const a = document.createElement('a', {});
    a.setAttribute("style", "display:none");
    a.setAttribute("href", url);
    a.setAttribute("download", file.name);
    document.body.appendChild(a);
    a.click();
    window.URL.revokeObjectURL(url);
}
function getAppConfiguration() {
    const filePicker = document.getElementById("filepicker");
    const fileList = filePicker.files;
    if ((fileList === null || fileList === void 0 ? void 0 : fileList.length) != 1)
        throw new Error("A single file must be selected");
    const inputFile = fileList[0];
    const alignInput = document.getElementById("align");
    const alignValue = alignInput.valueAsNumber;
    const truncateInput = document.getElementById("truncate");
    const truncateValue = truncateInput.checked;
    const orientationInput = document.querySelector('input[name="orientation"]:checked');
    const orientationValue = Orientation[orientationInput.value];
    const scoreFunctionInput = document.querySelector('input[name="scoringEngine"]:checked');
    const scoringFunctionValue = ScoreFunctions[scoreFunctionInput.value];
    const outputFormatInput = document.querySelector('input[name="outputFormat"]:checked');
    const outputFormatValue = OutputFormat[outputFormatInput.value];
    return {
        inputFile: inputFile,
        maxHeptadDisplacement: alignValue,
        truncate: truncateValue,
        orientation: orientationValue,
        scoreFunction: scoringFunctionValue,
        outputFormat: outputFormatValue
    };
}
function showCommandLine(args) {
    const cmdlineDiv = document.getElementById("fastscore-cmdline-container");
    cmdlineDiv === null || cmdlineDiv === void 0 ? void 0 : cmdlineDiv.removeAttribute("style");
    const cmdlineEl = document.getElementById("fastscore-cmdline");
    const cmdlineText = "fastscore " + args.join(" ");
    cmdlineEl.innerText = cmdlineText;
}
function runButtonClick() {
    return __awaiter(this, void 0, void 0, function* () {
        console.log("Button click");
        const appConfiguration = getAppConfiguration();
        const outputFiles = yield fastscoreCommandlineWrapper(appConfiguration);
        for (const file of outputFiles) {
            offerFileForDownload(file);
        }
    });
}
function initializeUI() {
    var _a;
    (_a = document.getElementById("run")) === null || _a === void 0 ? void 0 : _a.addEventListener("click", runButtonClick);
}
