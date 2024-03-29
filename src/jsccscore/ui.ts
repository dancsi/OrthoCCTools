

function offerFileForDownload(file: File): void {
    const url = window.URL.createObjectURL(file);

    const p = document.createElement('p');
    const a = document.createElement('a', {});
    p.append(a);

    a.setAttribute("href", url);
    a.setAttribute("download", file.name);

    const currentTime = new Date();
    a.innerText = `${file.name} (${currentTime.toLocaleString()})`;

    a.onclick = () => { 
        window.URL.revokeObjectURL(url); 
        a.setAttribute("class", "linkDisabled");
        a.onclick = () => {
            alert("You have already downloaded this file");
        }
    }

    const downloadAreaDiv = document.getElementById("downloadarea");
    downloadAreaDiv?.prepend(p);
    downloadAreaDiv?.removeAttribute("style");
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

    const outputFormatInput = document.querySelector('input[name="outputFormat"]:checked');
    const outputFormatValue = OutputFormat[(outputFormatInput as HTMLInputElement).value];

    return {
        inputFile: inputFile,
        maxHeptadDisplacement: alignValue,
        truncate: truncateValue,
        orientation: orientationValue,
        scoreFunction: scoringFunctionValue,
        outputFormat: outputFormatValue
    };
}

function showCommandLine(args: string[]): void
{
    const cmdlineDiv = document.getElementById("fastscore-cmdline-container");
    cmdlineDiv?.removeAttribute("style");
    
    const cmdlineEl = document.getElementById("fastscore-cmdline");
    const cmdlineText = "fastscore " + args.join(" ");

    cmdlineEl!.innerText = cmdlineText;
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