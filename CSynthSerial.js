// file CSynthSerial.js
// Load a series of xyz files from single drop, an handle them serially.
// Drag drop this text onto CSynth first,
// That will override the standard drag-drop hanling for the remainder of the session
// Then dragdrop the entire set of .xyz files
// Each one will be displayed as xyz positions without any dynamics, and saved as png image (to downloads)

let baseopenfiles;
const {openfiles, CSynth, log, V, GX, S, canvas, saveAs, W} = window;
baseopenfiles = baseopenfiles || openfiles;         // remember the original handler for use on each individual file
CSynth.applyXyzs = () => {};                        // don't bother with distance based model, minor optimization

/** new multi-file handler */
W.openfiles = async function openfilesserial(files) {
    for (const f of files) {                            // interactet the file list
        W.openfiles.pending = {};                       // for
        W.openfiles.dropped = {};                       // contents of dropped and other opened files
        await baseopenfiles([f])                        // standard process on single file
        await S.maestro('demoready')                    // wait for file etcs to be ready
        GX.getgui(/dynamicsrunning/).setValue(false)    // don't run dynamics
        GX.getgui(/matrix\/visible/).setValue(false)    // don't show matrix
	setBackgroundColor(1)
	
        await S.frame(2);

        GX.getgui(/positions/).press();         // go directly to positions
        GX.getgui(/autoscale/).press();         // and scale
        await S.frame(10);                      // small gap to settle (probably not needed, or only 1 or 2)
    }
}

W.openfiles.pending = {};
W.openfiles.dropped = {};  // contents of dropped and other opened files

