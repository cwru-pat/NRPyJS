/**
 * Spline class to help set initial data
 * To use: 
 *  - Include this file and Konva plots, <script src="https://unpkg.com/konva@8/konva.min.js"></script>
 *  - Add an element, <div id="chartContainer"></div>
 *  - Create a new spline instance, something like: var spl = new CSPL("chartContainer", 10);
 *  - Instead of makeInitials, call something like: spl.setSFID(NR, r, sf, psi, psi4, alpha, getScalarFieldID);
 **/ 

class CSPL
{
    constructor(chartContainer, n_pts) {

        // chart dimensions
        this.w = window.innerWidth/2;
        this.h = window.innerHeight/2;
        this.padding = 20;
        this.x_range = this.w - 2*this.padding;
        this.y_range = this.h - 2*this.padding;

        // Set up Konva stage
        this.stage = new Konva.Stage({
            container: chartContainer,
            width: this.w,
            height: this.h,
        });
        this.layer = new Konva.Layer();
        this.stage.add(this.layer);

        // draggable points
        this.xs = [], this.ys = [], this.draggable_points = [];
        this.n_pts = n_pts;
        for(var i=0; i<=n_pts; i++)
        {
            let y = this.padding + this.y_range/2 - this.y_range/2*Math.exp(-Math.pow(i/(n_pts*n_pts/40), 2));
            let x = this.padding + this.x_range*i/n_pts;
            this.xs[i] = x;
            this.ys[i] = y;
            this.draggable_points.push(this.draggablePoint(x, y, i, this.layer));
        }

        // interpolating function settings
        this.interp_n = 200;
        this.dx = 0.005*this.x_range/this.interp_n;
        this.ks = this.getNaturalKs(this.xs, this.ys);

        // add interpolating line, density line
        this.interpLine = new Konva.Line({
            points: this.getPtsList(),
            stroke: 'red',
            strokeWidth: 1,
            lineCap: 'round',
            lineJoin: 'round',
            tension: 1,
        });
        this.layer.add(this.interpLine);

        this.densityLine = new Konva.Line({
            points: this.getDensityPtsList(),
            stroke: '#074',
            strokeWidth: 1,
            lineCap: 'round',
            lineJoin: 'round',
            tension: 1,
        });
        this.layer.add(this.densityLine);
    }

    draggablePoint(x, y, pt_num, layer)
    {
        // create shape
        var box = new Konva.Circle({
            x: x,
            y: y,
            radius: 5,
            fill: '#C21',
            stroke: '#200',
            strokeWidth: 2,
            draggable: true,
        });
        layer.add(box);

        // constrain drag to y-direction only
        box.on('dragmove', () => {
            box.x(x);
            // box.y(Math.max(0, box.y));

            this.updatePaths();
        });

        // add cursor styling
        box.on('mouseover', function () {
            document.body.style.cursor = 'grab';
        });
        box.on('mouseout', function () {
            document.body.style.cursor = 'default';
        });

        return box;
    }

    updatePaths(box)
    {
        this.xs = [], this.ys = [];
        for(var i=0; i<=this.n_pts; i++) {
            let box = this.draggable_points[i];
            this.xs.push(box.attrs.x);
            this.ys.push(box.attrs.y);
        }
        this.ks = this.getNaturalKs(this.xs, this.ys);
        this.interpLine.setPoints(this.getPtsList());
        this.densityLine.setPoints(this.getDensityPtsList());
    }

    zerosMat(r, c) {
        var A = [];
        for(var i=0; i<r; i++) {
            A.push([]);
            for(var j=0; j<c; j++)
                A[i].push(0);
        }
        return A;
    }

    printMat(A) {
        for(var i=0; i<A.length; i++)
            console.log(A[i]);
    }

    swapRows(m, k, l) {
        var p = m[k];
        m[k] = m[l];
        m[l] = p;
    }

    solve(A, x) // in Matrix, out solutions
    {
        var m = A.length;
        for(var k=0; k<m; k++)  // column
        {
            // pivot for column
            var i_max = 0; var vali = Number.NEGATIVE_INFINITY;
            for(var i=k; i<m; i++) if(Math.abs(A[i][k])>vali) { i_max = i; vali = Math.abs(A[i][k]);}
            this.swapRows(A, k, i_max);
                        
            // for all rows below pivot
            for(var i=k+1; i<m; i++)
            {
                var cf = (A[i][k] / A[k][k]);
                for(var j=k; j<m+1; j++)  A[i][j] -= A[k][j] * cf;
            }
        }
        
        for(var i=m-1; i>=0; i--)   // rows = columns
        {
            var v = A[i][m] / A[i][i];
            x[i] = v;
            for(var j=i-1; j>=0; j--)   // rows
            {
                A[j][m] -= A[j][i] * v;
                A[j][i] = 0;
            }
        }
    }

    getNaturalKs(xs, ys) {
        var n = xs.length-1;
        var A = this.zerosMat(n+1, n+2);
        var ks = [];
        
        for(var i=1; i<n; i++)    // rows
        {
            A[i][i-1] = 1/(xs[i] - xs[i-1]);
            A[i][i  ] = 2 * (1/(xs[i] - xs[i-1]) + 1/(xs[i+1] - xs[i]));
            A[i][i+1] = 1/(xs[i+1] - xs[i]);
            A[i][n+1] = 3*   ( (ys[i]-ys[i-1])/ ((xs[i] - xs[i-1])*(xs[i] - xs[i-1])) 
                             +  (ys[i+1]-ys[i])/ ((xs[i+1] - xs[i])*(xs[i+1] - xs[i])) );
        }
        
        // zero first derivative at inner boundary?
        A[0][0  ] = 1.0;
        A[0][n+1] = 0.0;
        
        // Zero second-derivative at outer boundary
        A[n][n-1] = 1/(xs[n] - xs[n-1]);
        A[n][n  ] = 2/(xs[n] - xs[n-1]);
        A[n][n+1] = 3 * (ys[n] - ys[n-1]) / ((xs[n]-xs[n-1])*(xs[n]-xs[n-1]));
            
        this.solve(A, ks);

        return ks;
    }

    evalSpline(x, xs, ys, ks) {
        var i = 1;
        while(xs[i]<x) i++;
            
        var t = (x - xs[i-1]) / (xs[i] - xs[i-1]);
            
        var a =  ks[i-1]*(xs[i]-xs[i-1]) - (ys[i]-ys[i-1]);
        var b = -ks[i  ]*(xs[i]-xs[i-1]) + (ys[i]-ys[i-1]);
            
        var q = (1-t)*ys[i-1] + t*ys[i] + t*(1-t)*(a*(1-t)+b*t);
        return q;
    }

    evalAt(x) {
        return this.evalSpline(x, this.xs, this.ys, this.ks)
    }

    getPtsList() {
        let pts = [];
        for(let i=0; i<=2*this.interp_n; i+=2) {
            let x = this.padding + this.x_range*i/this.interp_n/2;
            pts[i] = x; // x value
            pts[i+1] = this.evalAt(x); // y value (same array)
        }
        return pts;
    }

    getDensityPtsList() {
        let pts = this.getPtsList();
        let derivs = [];
        for(let i=0; i<=2*this.interp_n; i+=2) {
            let x = this.padding + this.x_range*i/this.interp_n/2;
            let ptsp = this.evalAt(x-this.dx);
            derivs[i] = x;
            derivs[i+1] = this.h - this.padding - 5*Math.pow( (pts[i+1]-ptsp)/this.dx, 2 );
        }

        return derivs;
    }


    setSFID(NR, r, sf, psi, psi4, alpha,
        getScalarFieldID) {

        var rmax = 12.0;
        for(var i=0; i<NR; i++)
            r[i] = (i+0.5)/NR*rmax;

        // scalar field profile
        for(var i=0; i<NR; i++) {
            let x = this.padding + this.x_range*i/NR;
            sf[i] = this.evalAt(x) / this.h * 0.1;
        }

        getScalarFieldID(r, sf, psi);
    
        console.log('Generating Psi4 and Alpha! Test!');
        for(var i=0; i<NR; i++)
            psi4[i] = Math.pow(psi[i],4);
        for(var i=0; i<NR; i++)
            alpha[i] = Math.pow(psi[i],-2);

    }

}
