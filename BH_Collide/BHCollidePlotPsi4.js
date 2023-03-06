// Javascript file to plot simulation data

// Turns on/off debugging comments
var verbose = false;

// Various global variables
var runonce = false;
var fx, min, i0max, i1max, i2max, i0MIN, i1mid, i2mid;
var n, i, j, h, q, isRunning, firstRun, N_final, numStep, elem, barWidth;
var xdata=[], ydata=[], zdata=[], psi4XData=[], psi4YData=[];

// Updates the plot when different fields are selected
function changeUp(){
    if(runonce){
        fx = document.getElementById("selectFx").value;

        if(fx==24 || fx==25 || fx==26){
            updateMesh(zdata[fx][0], psi4XData[0], psi4YData[0]);
        } else {
            updateMesh(zdata[fx][0], xdata[0], ydata[0]);
        }
    }
};

function stopSim(){
    isRunning = false;
};

// Progress bar
function move(j) {
    if (j == 0) {
        elem = document.getElementById("progressDisplay");
        barWidth = 0;
    }
        if (barWidth >= 100) {
                clearInterval(id);
        } else {
            barWidth = barWidth + (100/(numStep+1));
            elem.style.width = barWidth + "%";
        }
};

// Main function that runs the simulation and records data into xdata,ydata,zdata
async function graphSlice(){
    if(n>N_final){
        isRunning = false;
    }

    if(isRunning){
        if(n%_getDim(3) == 0) {
            if(verbose) console.log("Generating output for plotting", n);

            i=0;
            for(let i2=min;i2<i2max;i2++){
                for(let i1=min;i1<i1max;i1++){
                    for(let i0=min;i0<i0max;i0++){
                        let idx = _getIDX3S(i0,i1,i2);
                        var xCart = [];
                        xCart[0] = _getCart(i0,i1,i2,1);
                        xCart[1] = _getCart(i0,i1,i2,2);
                        xdata[j][i] = xCart[0];
                        ydata[j][i] = xCart[1];
                        for(let k=0; k<24; k++){
                            zdata[k][j][i] = _getFxVal(_getIDX4ptS(k,idx));
                        }
                        if(verbose) console.log("Function: ", _getFxVal(_getIDX4ptS(fx,idx)));

                        i++;
                    }
                }
            }

            h=0;
            for(let r_ext_idx = 2;
                r_ext_idx<(_getNGHOSTS(1)-_getNGHOSTS(0))-2; r_ext_idx+=1) {

                let i0=r_ext_idx;
                for(let i1=min; i1<i1max; i1++) {
                    for(let i2=min; i2<i2max; i2++) {

                        if((i1 == (i1max-1)) && (i2 == (i2max-1))){
                            _runSpinWeightFx(n, r_ext_idx, i0, i1, i2, 1);
                        } else {
                            _runSpinWeightFx(n, r_ext_idx, i0, i1, i2, 0);
                        }

                        var xCart = [];
                        xCart[0] = _getCart(i0,i1,i2,1);
                        xCart[1] = _getCart(i0,i1,i2,2);
                        psi4XData[j][h] = xCart[0];
                        psi4YData[j][h] = xCart[1];

                        let idx3d = _getIDX3S(i0,i1,i2);
                        let psi4R = _getPsi4R(idx3d);
                        let psi4I = _getPsi4I(idx3d);
                        zdata[24][j][h] = Math.log10(Math.sqrt(psi4R**2)) + 6;
                        zdata[25][j][h] = Math.log10(Math.sqrt(psi4I**2)) + 6;
                        zdata[26][j][h] = Math.log10(Math.sqrt(psi4R**2 + psi4I**2)) + 6;

                        if(verbose) console.log("Psi4 data: ", xCart[0], xCart[1], psi4R, psi4I, zdata[26][j][f]);
                        h++;
                    }
                }
            }

            if(verbose) console.log("Done generating output for plotting", n);
            if(verbose) console.log("Updating Mesh");

            if(fx==24 || fx==25 || fx==26){
                updateMesh(zdata[fx][j], psi4XData[j], psi4YData[j]);
            } else {
                updateMesh(zdata[fx][j], xdata[j], ydata[j]);
            }

            await sleep(10);
            await move(j);
            j++;
        }

        if(verbose) console.log("Taking step", n);
        _stepForward();
        if(verbose) console.log("Done taking step", n);

        if(n==N_final){
            if(verbose) console.log("Generating output for plotting", n);
            i=0;
            for(let i2=min;i2<i2max;i2++){
                for(let i1=min;i1<i1max;i1++){
                    for(let i0=min;i0<i0max;i0++){
                        let idx = _getIDX3S(i0,i1,i2);
                        var xCart = [];
                        xCart[0] = _getCart(i0,i1,i2,1);
                        xCart[1] = _getCart(i0,i1,i2,2);
                        xdata[j][i] = xCart[0];
                        ydata[j][i] = xCart[1];
                        for(let funct=0; funct<24; funct++){
                            zdata[funct][j][i] = _getFxVal(_getIDX4ptS(funct,idx));
                        }
                        if(verbose) console.log(xCart[0], xCart[1], _getIDX4ptS(fx,idx), _getFxVal(_getIDX4ptS(fx,idx)));
                        i++;
                    }
                }
            }
            if(verbose) console.log("final ran", n);
            updateMesh(zdata[fx][j], xdata[j], ydata[j]);
            await sleep(10);
            j++;
            move();
        }
        n++;
        graphSlice();

    } else {
        document.getElementById("stopButton").disabled = true;
        document.getElementById("stopButton").style.cursor = "default";
        document.getElementById("nextButton").disabled = false;
        document.getElementById("nextButton").style.cursor = "pointer";
        document.getElementById("animateButton").disabled = false;
        document.getElementById("animateButton").style.cursor = "pointer";
        document.getElementById("progressBar").style.visibility = "hidden";
    }
};

//Moves plot forward one timestep in current function
function nextSlice(){
    if(q==j){
        q=0;
    }

    if(fx==24 || fx==25 || fx==26){
        updateMesh(zdata[fx][q], psi4XData[q], psi4YData[q]);
    } else {
        updateMesh(zdata[fx][q], xdata[q], ydata[q]);
    }
    q++;
}

//Animates through all timesteps for current function
async function animateSim(){
    q=0;
    for(var loop=0; loop<j; loop++){
        await nextSlice();
        await sleep(30);
    }
};

//Initializes variables and simulation
function runSim(){
    document.getElementById("stopButton").disabled = false;
    document.getElementById("stopButton").style.cursor = "pointer";
    document.getElementById("nextButton").disabled = true;
    document.getElementById("nextButton").style.cursor = "default";
    document.getElementById("animateButton").disabled = true;
    document.getElementById("animateButton").style.cursor = "default";
    document.getElementById("progressBar").style.visibility = "visible";

    _setDim(document.getElementById("res1").value,
            document.getElementById("res2").value,
            document.getElementById("res3").value);

    if(verbose) console.log("Running");
    _runsim();
    if(verbose) console.log("Done");

    N_final = _getNFinal();
    min = _getNGHOSTS(0);
    i0MIN=_getNGHOSTS(0); // In spherical, r=Delta r/2.
    i1mid= _getNGHOSTS(2)/2;
    i2mid= _getNGHOSTS(3)/2;
    i0max = _getNGHOSTS(1)-min;
    i1max = _getNGHOSTS(2)-min;
    i2max = _getNGHOSTS(3)-min;

    fx = document.getElementById("selectFx").value;
    isRunning = true;
    n=0;
    i=0;
    j=0;
    q=0;
    runonce=true;

    numStep = N_final/_getDim(3);

    for(let width=0; width<27; width++){
        zdata[width] = [];
        for(let length=0; length<(numStep+1); length++){
            xdata[length] = [];
            ydata[length] = [];
            psi4XData[length] = [];
            psi4YData[length] = [];
            zdata[width][length] = [];
        }
    }

    graphSlice();
};

function plotPsi4(){
    renderPsi4();
};

// Needed for plot to display in time
function sleep(ms){
    return new Promise(resolve => setTimeout(resolve, ms));
};
