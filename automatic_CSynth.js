var {log, sleep, openfiles, CSynth, savshowDirectoryPicker, Files, saveimage, G, S, monitorX, GX, V, distxyz, springs} = window;  // globals from javascript and CSynth


// rmse function on arrays of vectors
CSynth.rmsev = (a,b) => Math.sqrt(a.reduce((c,v,i) => c + distxyz(a[i], b[i]), 0))


/** process a single input txt file */
async function process(file, threshold = 0.1) {
    log('csyscript: realwork', file.name)
    await openfiles([await file.getFile()]);        // this will behave the same as dropping the single file
    await S.waitVal(() => CSynth.current.ready);    // wait till really started

    GX.getgui('ribbon/diameter').setValue(10);
    GX.getgui(/beddatasource/).setValue('rainbow');
    setBackgroundColor(1)
    G.stepsPerStep = 50       // 50 simulation steps per graphics display, less wasted graphics effort, more GPU for simulation
    G.springrate = 10       // bigger simulation steps for faster convergence
    G.springrate = 1        // standard simulation steps for better stability
    G.springpow = 0
    G.contactforce = 60
    G.pushapartpow = -3
    await sleep(300000) // divide by 1000 for wait time in seconds
    
    CSynth.showEigen(true);
    await S.frame()         // sleep till next frame (so eigen will have effect)
    CSynth.autoscale(0.1)
    await S.frame()        // sleep till next frame (so autoscale will have effect)
    G.scaleFactor = GX.getgui('modes/scaleFactor').getValue()*0.8;
    await S.frame()        // sleep till next frame (so zoom in/out will have effect)

    log('csyscript: realwork done, saving', file.name)
    CSynth.savexyz(file.name + '_3D.xyz')
    const sgui = V.gui.visible
    V.gui.visible = false;
    await saveimage(3000,2000, false, false, file.name + '.tga')
    V.gui.visible = sgui;
    log('csyscript: savingdone', file.name)
}


/** main test file */
async function test() {
    const dirhandle = Files.dirhandle || await Files.setDirectory();           // this will establish the Files local directory to work in
    const nts = dirhandle.entries();


    log('csyscript: starting');


    for await (const [name, file] of nts)
        if (name.endsWith('.txt')) await process(file);
    log('csyscript: all files done');
}

test()

// convergence loop
    // let last = springs.getpos();
    // let loop;
    // for (loop = 0; loop < 100; loop++) {     // loop till near enough
    //     await sleep(1000)
    //     const now = springs.getpos();
    //     const diff = CSynth.rmsev(last, now) ;
    //     log ('loop', loop, 'diff', diff);
    //     if (diff < threshold) break;
    //     last = now;
    // }
    // log('converged after loops', loop)
