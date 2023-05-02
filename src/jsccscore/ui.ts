

function offerFileForDownload(file: File): void {
    const url = window.URL.createObjectURL(file);

    const a = document.createElement('a', {});

    a.setAttribute("style", "display:none");
    a.setAttribute("href", url);
    a.setAttribute("download", file.name);

    document.body.appendChild(a);
    a.click();

    window.URL.revokeObjectURL(url);
}


function getAppConfiguration( ): FastscoreArguments {
    const filePicker = document.getElementById("filepicker");
    const fileList = (filePicker as HTMLInputElement).files;

    if(fileList?.length != 1) throw new Error("A single file must be selected");

    const inputFile = fileList[0];

    const alignInput = document.getElementById("align");
    const alignValue = (alignInput as HTMLInputElement).valueAsNumber;

    const truncateInput = document.getElementById("truncate");
    const truncateValue = (truncateInput as HTMLInputElement).checked;

    const orientationInput = document.querySelector('input[name="orientation"]:checked');
    const orientationValue = Orientation[(orientationInput as HTMLInputElement).value];

    const scoreFunctionInput = document.querySelector('input[name="scoringEngine"]:checked');
    const scoringFunctionValue = ScoreFunctions[(scoreFunctionInput as HTMLInputElement).value];

    return {
        inputFile: inputFile,
        maxHeptadDisplacement: alignValue,
        truncate: truncateValue,
        orientation: orientationValue,
        scoreFunction: scoringFunctionValue
    };
}

async function runButtonClick(): Promise<void>
{
    console.log("Button click");

    const appConfiguration = getAppConfiguration();

    const outputFiles = await fastscoreCommandlineWrapper(appConfiguration);

    for(const file of outputFiles) {
        offerFileForDownload(file);
    }
}

function initializeUI() {
    document.getElementById("run")?.addEventListener("click", runButtonClick);
}