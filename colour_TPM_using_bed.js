window.bedline2col = p => {
    let v = p[8]; //expect TPM in 9th column in bed file
    let threshold = 2;
    v = clamp(v, -threshold, threshold);
    v = v/threshold;  //map to [-1,1]
    let r,g,b;
    if (v < 0) {                  // blue to white to red diverging colourmap
          b = 1; r = g = 1+v;
	  //r = 1; g = b = 1; //white
    } else {
          r = 1; g = b = 1-v;
    }
    return {r, g, b};
}
