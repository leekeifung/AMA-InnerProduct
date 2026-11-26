/**
 * AMA1751 Linear Algebra — Topic 6: Inner Product
 * Module 9: Applied Orthogonality — Numerical Stability (CGS vs MGS vs Householder),
 *           Cosine Similarity in AI, and Orthonormal Signal Denoising (DCT & Parseval)
 */
(function () {
  function fmtSci(x) {
    if (!Number.isFinite(x)) return 'NaN';
    if (x === 0) return '0.00e+0';
    if (x >= 1e-2 && x < 1e4) return x.toFixed(4);
    return x.toExponential(2);
  }

  function clean(x, d = 3) {
    if (!Number.isFinite(x)) return 'NaN';
    if (Math.abs(x) < 1e-12) return '0';
    return String(Number(x.toFixed(d)));
  }

  function dot(a, b) {
    let s = 0;
    for (let i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }

  function norm(a) {
    return Math.sqrt(dot(a, a));
  }

  // =========================================================================
  // LAB 1 MATH: LÄUCHLI MATRIX & REAL IEEE-754 CGS / MGS / HOUSEHOLDER
  // =========================================================================
  function lauchliCols(eps) {
    return [
      [1, eps, 0, 0],
      [1, 0, eps, 0],
      [1, 0, 0, eps]
    ];
  }

  function classicalGS(cols) {
    const k = cols.length;
    const m = cols[0].length;
    const Q = [];
    for (let j = 0; j < k; j++) {
      const v = cols[j].slice();
      // CGS computes all inner products against the ORIGINAL column cols[j]
      const r = [];
      for (let i = 0; i < j; i++) {
        r.push(dot(Q[i], cols[j]));
      }
      for (let i = 0; i < j; i++) {
        for (let row = 0; row < m; row++) v[row] -= r[i] * Q[i][row];
      }
      const nv = norm(v);
      Q.push(nv > 0 ? v.map(x => x / nv) : v);
    }
    return Q;
  }

  function modifiedGS(cols) {
    const k = cols.length;
    const m = cols[0].length;
    const Q = [];
    for (let j = 0; j < k; j++) {
      const v = cols[j].slice();
      // MGS updates v sequentially before taking the next inner product
      for (let i = 0; i < j; i++) {
        const rij = dot(Q[i], v);
        for (let row = 0; row < m; row++) v[row] -= rij * Q[i][row];
      }
      const nv = norm(v);
      Q.push(nv > 0 ? v.map(x => x / nv) : v);
    }
    return Q;
  }

  function householderQR(cols) {
    const n = cols.length;    // 3
    const m = cols[0].length; // 4
    // Working copy A[row][col]
    const A = Array.from({ length: m }, (_, r) => cols.map(c => c[r]));
    // Q starts as m x m identity
    const Qfull = Array.from({ length: m }, (_, r) =>
      Array.from({ length: m }, (__, c) => (r === c ? 1 : 0))
    );

    for (let k = 0; k < n; k++) {
      const x = [];
      for (let r = k; r < m; r++) x.push(A[r][k]);
      const nx = norm(x);
      if (nx === 0) continue;
      const sign = x[0] >= 0 ? 1 : -1;
      const u = x.slice();
      u[0] += sign * nx;
      const uNormSq = dot(u, u);
      if (uNormSq === 0) continue;
      const beta = 2 / uNormSq;

      // Apply H_k to A[k:m, k:n] from left
      for (let c = k; c < n; c++) {
        let s = 0;
        for (let r = 0; r < m - k; r++) s += u[r] * A[k + r][c];
        for (let r = 0; r < m - k; r++) A[k + r][c] -= beta * s * u[r];
      }

      // Accumulate Qfull = Qfull * H_k (apply H_k to columns k:m from right)
      for (let r = 0; r < m; r++) {
        let s = 0;
        for (let c = 0; c < m - k; c++) s += Qfull[r][k + c] * u[c];
        for (let c = 0; c < m - k; c++) Qfull[r][k + c] -= beta * s * u[c];
      }
    }

    // Extract first n columns of Qfull
    const Q = [];
    for (let c = 0; c < n; c++) {
      const col = [];
      for (let r = 0; r < m; r++) col.push(Qfull[r][c]);
      Q.push(col);
    }
    return Q;
  }

  function gramError(Q) {
    const k = Q.length;
    const E = Array.from({ length: k }, () => Array(k).fill(0));
    let frobSq = 0;
    for (let i = 0; i < k; i++) {
      for (let j = 0; j < k; j++) {
        const val = dot(Q[i], Q[j]) - (i === j ? 1 : 0);
        E[i][j] = Math.abs(val);
        frobSq += val * val;
      }
    }
    return { E, frob: Math.sqrt(frobSq) };
  }

  // =========================================================================
  // LAB 2 DATA: SEMANTIC VECTOR SEARCH DOCUMENTS (R^4)
  // =========================================================================
  const DOCS = [
    { id: 'A', name: 'Doc A · Orthogonal Projections & QR',   vec: [9, 1, 3, 2], color: '#0F766E' },
    { id: 'B', name: 'Doc B · Least Squares & Regression',    vec: [7, 2, 6, 8], color: '#B45309' },
    { id: 'C', name: 'Doc C · Deep Learning Embeddings',      vec: [6, 4, 9, 7], color: '#1D4ED8' },
    { id: 'D', name: 'Doc D · Multivariable Optimization',    vec: [3, 9, 4, 2], color: '#7C3AED' },
    { id: 'E', name: 'Doc E · Probability & Bayesian Stats',  vec: [2, 5, 4, 9], color: '#BE123C' }
  ];

  const QUERY_PRESETS = {
    la:   [9, 1, 2, 1],
    ds:   [7, 2, 6, 8],
    ml:   [5, 4, 9, 7],
    calc: [2, 9, 3, 2]
  };

  // =========================================================================
  // LAB 3 MATH: ORTHONORMAL DCT-II BASIS (N = 64) & SIGNALS
  // =========================================================================
  const SIG_N = 64;
  const DCT_BASIS = (() => {
    const basis = [];
    for (let k = 0; k < SIG_N; k++) {
      const alpha = k === 0 ? Math.sqrt(1 / SIG_N) : Math.sqrt(2 / SIG_N);
      const v = new Array(SIG_N);
      for (let n = 0; n < SIG_N; n++) {
        v[n] = alpha * Math.cos((Math.PI * k * (n + 0.5)) / SIG_N);
      }
      basis.push(v);
    }
    return basis;
  })();

  function makeCleanSignal(type) {
    const y = new Array(SIG_N);
    for (let n = 0; n < SIG_N; n++) {
      const t = (n + 0.5) / SIG_N;
      if (type === 'multitone') {
        y[n] = 1.0 * Math.cos(Math.PI * t) + 0.65 * Math.cos(3 * Math.PI * t) - 0.35 * Math.cos(5 * Math.PI * t);
      } else if (type === 'ecg') {
        const pWave = 0.25 * Math.exp(-Math.pow((t - 0.22) / 0.05, 2));
        const qDip  = -0.30 * Math.exp(-Math.pow((t - 0.38) / 0.02, 2));
        const rPeak = 1.65 * Math.exp(-Math.pow((t - 0.42) / 0.022, 2));
        const sDip  = -0.40 * Math.exp(-Math.pow((t - 0.46) / 0.02, 2));
        const tWave = 0.45 * Math.exp(-Math.pow((t - 0.68) / 0.075, 2));
        y[n] = pWave + qDip + rPeak + sDip + tWave - 0.15;
      } else if (type === 'step') {
        y[n] = t < 0.5 ? 1.0 : -1.0;
      } else {
        y[n] = 2.6 * Math.pow(t - 0.5, 2) - 0.4 + 0.45 * Math.sin(6 * Math.PI * t);
      }
    }
    return y;
  }

  function makeGaussianNoise(seed) {
    // Deterministic Box-Muller via Mulberry32 PRNG
    let a = seed >>> 0;
    function rand() {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    }
    const z = new Array(SIG_N);
    for (let i = 0; i < SIG_N; i += 2) {
      const u1 = Math.max(1e-9, rand());
      const u2 = rand();
      const r = Math.sqrt(-2 * Math.log(u1));
      z[i] = r * Math.cos(2 * Math.PI * u2);
      if (i + 1 < SIG_N) z[i + 1] = r * Math.sin(2 * Math.PI * u2);
    }
    return z;
  }

  // =========================================================================
  // CAPSTONE CONCEPTUAL CHECKS
  // =========================================================================
  const CAPSTONE_QS = [
    {
      q: 'In exact arithmetic, Classical Gram–Schmidt (CGS) and Modified Gram–Schmidt (MGS) produce identical $Q$ matrices. Why does CGS fail catastrophically on the Läuchli matrix when $\\varepsilon \\approx 10^{-8}$ in IEEE-754 double precision?',
      opts: [
        'Because CGS divides by zero whenever $\\varepsilon < 10^{-4}$.',
        'Because CGS projects the original column $\\mathbf{a}_j$ against all previous $\\mathbf{q}_i$ simultaneously, so roundoff error from $\\mathbf{q}_1$ is never removed before taking $\\mathbf{q}_2\\cdot\\mathbf{a}_3$.',
        'Because MGS uses Householder reflections internally.',
        'Because the Läuchli matrix has linearly dependent columns for all $\\varepsilon < 1$.'
      ],
      ans: 1,
      exp: 'In CGS, $\\mathbf{q}_2\\cdot\\mathbf{a}_3$ relies on the severely cancelled first entry of $\\mathbf{q}_2$ (which rounds to $0$ at $\\varepsilon=10^{-8}$). MGS first subtracts the $\\mathbf{q}_1$-component from $\\mathbf{a}_3$, extracting the exact component along $\\mathbf{q}_2$ from the uncancelled entries!'
    },
    {
      q: 'In a vector search engine, Document B is concatenated with a duplicate copy of itself ($\\mathbf{d}_B \\mapsto 2\\mathbf{d}_B$). How do its Cosine Similarity and Euclidean Distance to a fixed query $\\mathbf{q}$ behave?',
      opts: [
        'Both Cosine Similarity and Euclidean Distance remain unchanged.',
        'Cosine Similarity doubles, while Euclidean Distance stays unchanged.',
        'Cosine Similarity stays identical (scale-invariant), whereas Euclidean Distance changes and can distort search rankings.',
        'Cosine Similarity drops to $0$, while Euclidean Distance halves.'
      ],
      ans: 2,
      exp: 'Since $\\frac{\\mathbf{q}\\cdot(s\\mathbf{d})}{\\|\\mathbf{q}\\|\\|s\\mathbf{d}\\|} = \\frac{s(\\mathbf{q}\\cdot\\mathbf{d})}{s\\|\\mathbf{q}\\|\\|\\mathbf{d}\\|} = \\cos\\theta$ for any $s > 0$, cosine similarity depends purely on direction (semantic topic mix), not document length.'
    },
    {
      q: 'Let $\\{\\mathbf{q}_0, \\dots, \\mathbf{q}_{N-1}\\}$ be an orthonormal basis of $\\mathbb{R}^N$ (such as DCT-II) and $c_k = \\mathbf{y}\\cdot\\mathbf{q}_k$. If we project $\\mathbf{y}$ onto the first $K$ basis vectors ($\\hat{\\mathbf{y}}_K = \\sum_{k=0}^{K-1} c_k \\mathbf{q}_k$), what is the exact squared compression error $\\|\\mathbf{y} - \\hat{\\mathbf{y}}_K\\|^2$?',
      opts: [
        '$\\sum_{k=0}^{K-1} c_k^2$',
        '$\\sum_{k=K}^{N-1} c_k^2$ (the energy of the discarded coordinates, by Parseval’s Identity)',
        '$\\left(\\sum_{k=K}^{N-1} |c_k|\\right)^2$',
        '$0$, whenever $K \\ge N/2$'
      ],
      ans: 1,
      exp: 'Because $\\mathbf{y} - \\hat{\\mathbf{y}}_K = \\sum_{k=K}^{N-1} c_k \\mathbf{q}_k$ and the $\\mathbf{q}_k$ are orthonormal, Pythagoras/Parseval gives $\\|\\mathbf{y} - \\hat{\\mathbf{y}}_K\\|^2 = \\sum_{k=K}^{N-1} c_k^2$.'
    }
  ];

  // =========================================================================
  // MODULE 9 REGISTRATION
  // =========================================================================
  AMA1751.registerModule(9, {
    render(meta) {
      const slideLabel = (meta && meta.slides) ? meta.slides : 'Extensions · Data Science, Numerical Computing & Signals';
      return `
        <style>
          .m9-controls {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            margin-bottom: 12px;
            align-items: center;
          }
          .m9-chip {
            padding: 6px 12px;
            border-radius: 999px;
            border: 1px solid var(--border);
            background: #F8FAFC;
            font-size: 0.82rem;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.15s;
          }
          .m9-chip:hover, .m9-chip.active {
            background: var(--polyu-red-soft);
            border-color: #E8B4C4;
            color: var(--polyu-red-dark);
          }
          .m9-btn {
            padding: 7px 13px;
            border-radius: 8px;
            border: 1px solid #CBD5E1;
            background: #FFFFFF;
            font-weight: 700;
            font-size: 0.84rem;
            cursor: pointer;
          }
          .m9-btn:hover { background: #F1F5F9; }
          .m9-panel {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 14px 16px;
            font-size: 0.89rem;
            line-height: 1.55;
          }
          .m9-canvas-box {
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 8px;
          }
          .m9-canvas {
            display: block;
            width: 100%;
          }
          .m9-slider-row {
            display: flex;
            align-items: center;
            gap: 10px;
            margin: 6px 0;
            font-size: 0.85rem;
          }
          .m9-slider-row label {
            min-width: 130px;
            font-weight: 700;
            color: #334155;
          }
          .m9-slider-row input[type="range"] {
            flex: 1;
            accent-color: #0F766E;
          }
          .m9-slider-val {
            min-width: 72px;
            text-align: right;
            font-weight: 800;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
            color: #0F172A;
          }
          .m9-stat-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(130px, 1fr));
            gap: 8px;
            margin: 10px 0;
          }
          .m9-stat-card {
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 8px 10px;
            text-align: center;
          }
          .m9-stat-lbl {
            font-size: 0.69rem;
            font-weight: 700;
            text-transform: uppercase;
            color: var(--text-muted);
          }
          .m9-stat-val {
            font-size: 0.96rem;
            font-weight: 800;
            margin-top: 2px;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
            color: #0F172A;
          }
          .m9-table {
            width: 100%;
            border-collapse: collapse;
            font-size: 0.82rem;
            background: #FFFFFF;
            border-radius: 8px;
            overflow: hidden;
            border: 1px solid var(--border);
          }
          .m9-table th, .m9-table td {
            padding: 7px 8px;
            border-bottom: 1px solid #E2E8F0;
            text-align: left;
          }
          .m9-table th {
            background: #F1F5F9;
            font-weight: 700;
            color: #334155;
          }
          .m9-pill {
            display: inline-block;
            padding: 2px 7px;
            border-radius: 999px;
            font-size: 0.74rem;
            font-weight: 800;
          }
          .m9-pill.top { background: #DCFCE7; color: #166534; }
          .m9-pill.warn { background: #FEF3C7; color: #92400E; }
          .m9-pill.bad { background: #FEE2E2; color: #991B1B; }
          .m9-q-card {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 12px 14px;
            margin-bottom: 10px;
          }
          .m9-q-opt {
            display: block;
            width: 100%;
            text-align: left;
            padding: 8px 11px;
            margin-top: 6px;
            border: 1px solid #CBD5E1;
            border-radius: 8px;
            background: #FFFFFF;
            cursor: pointer;
            font-size: 0.86rem;
          }
          .m9-q-opt:hover { background: #F1F5F9; }
          .m9-q-opt.correct { border-color: #16A34A; background: #DCFCE7; color: #166534; font-weight: 700; }
          .m9-q-opt.wrong { border-color: #DC2626; background: #FEE2E2; color: #991B1B; }
        </style>

        <!-- HEADER -->
        <div class="module-header">
          <span class="slide-ref">Part 9 · ${slideLabel}</span>
          <h2>9. Applied Orthogonality: Numerical Stability, AI Search &amp; Signal Denoising</h2>
          <p>See how <strong>Topic 6 (Inner Product, Orthogonality, QR &amp; Projections)</strong> powers modern scientific computing, AI vector search, and audio/image compression—complete with live IEEE-754 experiments.</p>
        </div>

        <!-- =============================================================== -->
        <!-- LAB 1: NUMERICAL STABILITY (CGS vs MGS vs HOUSEHOLDER)          -->
        <!-- =============================================================== -->
        <div class="card">
          <h3>Interactive Lab 1: Classical GS vs. Modified GS vs. Householder QR (Live IEEE-754 Experiment)</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            On paper, Classical Gram–Schmidt (CGS) and Modified Gram–Schmidt (MGS) are identical. In real 64-bit floating-point arithmetic ($\varepsilon_{\\text{mach}} \\approx 2.22\\times 10^{-16}$), test them on the famous $4\\times 3$ <strong>Läuchli matrix</strong>:
          </p>
          <div style="overflow-x:auto; margin-bottom:10px;">
            $$A(\\varepsilon) = \\begin{bmatrix} 1 & 1 & 1 \\\\ \\varepsilon & 0 & 0 \\\\ 0 & \\varepsilon & 0 \\\\ 0 & 0 & \\varepsilon \\end{bmatrix}, \\qquad \\kappa(A) = \\frac{\\sqrt{3 + \\varepsilon^2}}{\\varepsilon}, \\qquad \\kappa(A^T A) = \\kappa(A)^2$$
          </div>

          <div class="m9-panel" style="margin-bottom:12px;">
            <div class="m9-slider-row">
              <label for="m9EpsSlider">Perturbation $\\varepsilon = 10^{-p}$:</label>
              <input type="range" id="m9EpsSlider" min="1" max="8" step="0.25" value="6" />
              <span class="m9-slider-val" id="m9EpsVal">1.00e-6</span>
            </div>
            <div class="m9-controls" style="margin:6px 0 0;">
              <strong style="font-size:0.8rem;">Jump to regime:</strong>
              <button class="m9-chip" data-eps="2">Well-Conditioned ($10^{-2}$)</button>
              <button class="m9-chip active" data-eps="6">Moderate Loss ($10^{-6}$)</button>
              <button class="m9-chip" data-eps="7.5">Critical Threshold ($10^{-7.5}$)</button>
              <button class="m9-chip" data-eps="8">CGS Collapse ($10^{-8} \\approx \\sqrt{\\varepsilon_{\\text{mach}}}$!)</button>
            </div>
          </div>

          <div class="grid-2" style="align-items:start;">
            <div class="m9-canvas-box">
              <canvas id="m9StabCanvas" class="m9-canvas" style="height:265px;"></canvas>
            </div>
            <div class="m9-panel">
              <div style="font-weight:700; margin-bottom:6px;">Real-Time IEEE-754 Orthogonality Diagnostics:</div>
              <div class="m9-stat-grid" id="m9StabStats"></div>
              <div id="m9StabVerdict" style="font-size:0.85rem; margin-top:8px;"></div>
            </div>
          </div>
        </div>

        <!-- =============================================================== -->
        <!-- LAB 2: COSINE SIMILARITY & AI SEMANTIC VECTOR SEARCH            -->
        <!-- =============================================================== -->
        <div class="card">
          <h3>Interactive Lab 2: Cosine Similarity vs. Euclidean Distance (AI Embeddings &amp; Search)</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            Modern AI retrieval (RAG, recommender systems, LLMs) ranks documents by <strong>Cosine Similarity</strong>—literally the inner product of unit-normalized vectors:
            $\\cos\\theta = \\frac{\\mathbf{q}\\cdot\\mathbf{d}}{\\|\\mathbf{q}\\|\\|\\mathbf{d}\\|} = \\hat{\\mathbf{q}}\\cdot\\hat{\\mathbf{d}}$.
            Drag the <strong>Doc B Length Multiplier $s$</strong> below to watch why Euclidean distance fails when document lengths vary!
          </p>

          <div class="m9-controls">
            <strong style="font-size:0.82rem;">Query Presets:</strong>
            <button class="m9-chip" data-qp="la">Pure Linear Algebra</button>
            <button class="m9-chip active" data-qp="ds">Data Science / Regression</button>
            <button class="m9-chip" data-qp="ml">AI / Deep Learning</button>
            <button class="m9-chip" data-qp="calc">Calculus &amp; Optimization</button>
          </div>

          <div class="grid-2" style="align-items:start;">
            <div class="m9-panel">
              <div style="font-weight:700; margin-bottom:6px;">Query Vector $\\mathbf{q} \\in \\mathbb{R}^4$ (Topic Weights):</div>
              <div class="m9-slider-row"><label>Linear Algebra ($q_1$):</label><input type="range" id="m9Q0" min="0" max="10" step="1" value="7" /><span class="m9-slider-val" id="m9Q0v">7</span></div>
              <div class="m9-slider-row"><label>Calculus ($q_2$):</label><input type="range" id="m9Q1" min="0" max="10" step="1" value="2" /><span class="m9-slider-val" id="m9Q1v">2</span></div>
              <div class="m9-slider-row"><label>Python / ML ($q_3$):</label><input type="range" id="m9Q2" min="0" max="10" step="1" value="6" /><span class="m9-slider-val" id="m9Q2v">6</span></div>
              <div class="m9-slider-row"><label>Statistics ($q_4$):</label><input type="range" id="m9Q3" min="0" max="10" step="1" value="8" /><span class="m9-slider-val" id="m9Q3v">8</span></div>

              <hr style="border:none; border-top:1px solid #CBD5E1; margin:10px 0;" />
              <div class="m9-slider-row">
                <label style="color:#B45309;">Doc B Length Scale ($s$):</label>
                <input type="range" id="m9DocScale" min="0.5" max="4.0" step="0.1" value="1.0" />
                <span class="m9-slider-val" id="m9DocScaleV" style="color:#B45309;">1.0×</span>
              </div>
              <div style="font-size:0.8rem; color:var(--text-muted); margin-top:4px;">
                Simulates making Doc B $s\\times$ longer (repeating the same topic proportions: $\\mathbf{d}_B \\mapsto s\\mathbf{d}_B$).
              </div>

              <div class="m9-canvas-box" style="margin-top:10px;">
                <canvas id="m9CosCanvas" class="m9-canvas" style="height:220px;"></canvas>
              </div>
            </div>

            <div>
              <table class="m9-table">
                <thead>
                  <tr>
                    <th>Document</th>
                    <th>Angle $\\theta$</th>
                    <th>$\\cos\\theta$ (Rank)</th>
                    <th>$\\|\\mathbf{q}-\\mathbf{d}\\|$ (Rank)</th>
                  </tr>
                </thead>
                <tbody id="m9CosTableBody"></tbody>
              </table>
              <div class="m9-panel" id="m9CosInsight" style="margin-top:10px; font-size:0.85rem;"></div>
            </div>
          </div>
        </div>

        <!-- =============================================================== -->
        <!-- LAB 3: ORTHONORMAL DCT COMPRESSION & SIGNAL DENOISING           -->
        <!-- =============================================================== -->
        <div class="card">
          <h3>Interactive Lab 3: Orthonormal Projection Denoising &amp; Parseval’s Identity (DCT-II Basis)</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            The $64$ Discrete Cosine Transform vectors $\\{\\mathbf{q}_0, \\dots, \\mathbf{q}_{63}\\}$ form an <strong>orthonormal basis</strong> of $\\mathbb{R}^{64}$ (the exact engine behind JPEG and MP3!). Projecting a noisy signal $\\mathbf{y}$ onto the low-frequency subspace $W_K = \\operatorname{span}\\{\\mathbf{q}_0,\\dots,\\mathbf{q}_{K-1}\\}$ via $\\hat{\\mathbf{y}}_K = \\sum_{k=0}^{K-1}(\\mathbf{y}\\cdot\\mathbf{q}_k)\\mathbf{q}_k$ removes high-frequency noise!
          </p>

          <div class="m9-controls">
            <strong style="font-size:0.82rem;">Signal:</strong>
            <button class="m9-chip active" data-sig="multitone">Smooth Multi-Tone</button>
            <button class="m9-chip" data-sig="ecg">ECG Heartbeat Pulse</button>
            <button class="m9-chip" data-sig="quad">Quadratic + Ripple</button>
            <button class="m9-chip" data-sig="step">Step Edge (Gibbs Test)</button>
            <button class="m9-btn" id="m9ReseedBtn" style="margin-left:auto;">↻ Resample Noise</button>
          </div>

          <div class="m9-panel" style="margin-bottom:12px;">
            <div class="m9-slider-row">
              <label for="m9KSlider">Kept Basis Modes ($K / 64$):</label>
              <input type="range" id="m9KSlider" min="1" max="64" step="1" value="8" />
              <span class="m9-slider-val" id="m9KVal">K = 8 (12.5%)</span>
            </div>
            <div class="m9-slider-row">
              <label for="m9NoiseSlider">Noise Level ($\\sigma$):</label>
              <input type="range" id="m9NoiseSlider" min="0" max="0.8" step="0.05" value="0.30" />
              <span class="m9-slider-val" id="m9NoiseVal">0.30</span>
            </div>
          </div>

          <div class="m9-canvas-box" style="margin-bottom:10px;">
            <canvas id="m9SigCanvas" class="m9-canvas" style="height:290px;"></canvas>
          </div>

          <div class="m9-stat-grid" id="m9SigStats"></div>

          <details class="code-details">
            <summary>Python (NumPy / SciPy) Companion Code — Stability, Cosine Search &amp; DCT Denoising</summary>
            <pre><code>import numpy as np
from scipy.fft import dct, idct

# 1) Numerical Stability on 4x3 Lauchli Matrix (eps = 1e-8)
eps = 1e-8
A = np.array([[1., 1., 1.], [eps, 0., 0.], [0., eps, 0.], [0., 0., eps]])
Q_hh, R_hh = np.linalg.qr(A, mode='reduced')  # Householder QR
print("Householder ||Q^T Q - I||_F =", np.linalg.norm(Q_hh.T @ Q_hh - np.eye(3)))
print("cond(A) =", f"{np.linalg.cond(A):.2e}", "| cond(A^T A) =", f"{np.linalg.cond(A.T @ A):.2e}")

# 2) Scale-Invariance of Cosine Similarity
q = np.array([7., 2., 6., 8.])
d = np.array([7., 2., 6., 8.])
cos_sim = lambda u, v: np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
print("cos(q, d) =", cos_sim(q, d), "| cos(q, 3*d) =", cos_sim(q, 3 * d))

# 3) Orthonormal DCT-II Projection Denoising + Parseval Check
N, K = 64, 8
t = (np.arange(N) + 0.5) / N
y_clean = np.cos(np.pi * t) + 0.65 * np.cos(3 * np.pi * t) - 0.35 * np.cos(5 * np.pi * t)
y_noisy = y_clean + 0.30 * np.random.default_rng(42).standard_normal(N)
c = dct(y_noisy, norm='ortho')                # Exact orthonormal coordinates c_k = y . q_k
print("Parseval check: ||y||^2 =", np.dot(y_noisy, y_noisy), "| sum(c^2) =", np.dot(c, c))
c_trunc = np.zeros_like(c); c_trunc[:K] = c[:K]
y_hat = idct(c_trunc, norm='ortho')           # Orthogonal projection onto W_K
print("Noisy error =", np.linalg.norm(y_noisy - y_clean), "| Denoised error =", np.linalg.norm(y_hat - y_clean))</code></pre>
          </details>
        </div>

        <!-- =============================================================== -->
        <!-- SECTION D: CAPSTONE CONCEPTUAL CHECK                            -->
        <!-- =============================================================== -->
        <div class="card">
          <h3>Applied Orthogonality Checkpoint (3 Conceptual Questions)</h3>
          <div id="m9CapstoneQuiz"></div>
        </div>
      `;
    },

    mount(container, { renderMath }) {
      // =====================================================================
      // LAB 1: NUMERICAL STABILITY (CGS vs MGS vs HOUSEHOLDER)
      // =====================================================================
      const epsSlider = container.querySelector('#m9EpsSlider');
      const epsValEl = container.querySelector('#m9EpsVal');
      const stabCanvas = container.querySelector('#m9StabCanvas');
      const stabCtx = stabCanvas.getContext('2d');

      function updateStabilityLab() {
        const p = parseFloat(epsSlider.value);
        const eps = Math.pow(10, -p);
        epsValEl.textContent = eps.toExponential(2);

        const cols = lauchliCols(eps);
        const qCGS = classicalGS(cols);
        const qMGS = modifiedGS(cols);
        const qHH  = householderQR(cols);

        const eCGS = gramError(qCGS);
        const eMGS = gramError(qMGS);
        const eHH  = gramError(qHH);

        const condA = Math.sqrt(3 + eps * eps) / eps;
        const condATA = condA * condA;

        // Stats cards
        container.querySelector('#m9StabStats').innerHTML = `
          <div class="m9-stat-card">
            <div class="m9-stat-lbl">CGS $\\|Q^TQ-I\\|_F$</div>
            <div class="m9-stat-val" style="color:${eCGS.frob > 1e-3 ? '#DC2626' : '#B45309'};">${fmtSci(eCGS.frob)}</div>
          </div>
          <div class="m9-stat-card">
            <div class="m9-stat-lbl">MGS $\\|Q^TQ-I\\|_F$</div>
            <div class="m9-stat-val" style="color:#1D4ED8;">${fmtSci(eMGS.frob)}</div>
          </div>
          <div class="m9-stat-card">
            <div class="m9-stat-lbl">Householder $\\|Q^TQ-I\\|_F$</div>
            <div class="m9-stat-val" style="color:#0F766E;">${fmtSci(eHH.frob)}</div>
          </div>
          <div class="m9-stat-card">
            <div class="m9-stat-lbl">$\\kappa(A)$ / $\\kappa(A^T A)$</div>
            <div class="m9-stat-val" style="font-size:0.82rem;">${fmtSci(condA)} / ${fmtSci(condATA)}</div>
          </div>
        `;
        renderMath(container.querySelector('#m9StabStats'));

        const verdictEl = container.querySelector('#m9StabVerdict');
        if (eCGS.frob > 0.1) {
          verdictEl.innerHTML = `<span class="m9-pill bad">CGS CATASTROPHIC COLLAPSE!</span> Because $1 + \\varepsilon^2$ rounds to $1.0$ in 64-bit float, CGS produces $\\mathbf{q}_2 \\cdot \\mathbf{q}_3 \\approx ${clean(dot(qCGS[1], qCGS[2]), 3)}$ instead of $0$! Meanwhile <strong>Householder QR</strong> stays at $\\approx 10^{-16}$, and $\\kappa(A^T A) \\approx ${fmtSci(condATA)}$ exceeds $1/\\varepsilon_{\\text{mach}}$ (normal equations lose all 16 digits!).`;
        } else if (eCGS.frob > 1e-6) {
          verdictEl.innerHTML = `<span class="m9-pill warn">Orthogonality Drifting</span> CGS error scales as $O(\\kappa(A)^2\\varepsilon_{\\text{mach}})$ while MGS scales as $O(\\kappa(A)\\varepsilon_{\\text{mach}})$. Notice MGS is already <strong>${fmtSci(eCGS.frob / Math.max(1e-17, eMGS.frob))}×</strong> more accurate than CGS!`;
        } else {
          verdictEl.innerHTML = `<span class="m9-pill top">Well-Conditioned Regime</span> All three algorithms preserve orthogonality well here. Push $\\varepsilon$ toward $10^{-7.5}$–$10^{-8}$ to watch Classical GS shatter!`;
        }
        renderMath(verdictEl);

        drawStabilityCanvas(eCGS, eMGS, eHH);
      }

      function drawStabilityCanvas(eCGS, eMGS, eHH) {
        if (!stabCanvas.isConnected) return;
        const dpr = window.devicePixelRatio || 1;
        const rect = stabCanvas.getBoundingClientRect();
        const w = rect.width || 440, h = rect.height || 265;
        stabCanvas.width = Math.round(w * dpr);
        stabCanvas.height = Math.round(h * dpr);
        stabCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
        stabCtx.clearRect(0, 0, w, h);

        // Left half: Log10 bar chart of ||Q^T Q - I||_F from 1e-16 to 1e0
        const leftW = w * 0.56;
        const padL = 42, padT = 26, padB = 34, chartH = h - padT - padB;
        stabCtx.fillStyle = '#0F172A';
        stabCtx.font = 'bold 11.5px sans-serif';
        stabCtx.textAlign = 'left';
        stabCtx.fillText('Log₁₀ Orthogonality Error ||QᵀQ − I||_F', 10, 16);

        // Gridlines at 10^0, 10^-4, 10^-8, 10^-12, 10^-16
        stabCtx.font = '10px monospace';
        stabCtx.textAlign = 'right';
        [0, -4, -8, -12, -16].forEach(exp => {
          const frac = (0 - exp) / 16; // 0 at top, 1 at bottom
          const y = padT + frac * chartH;
          stabCtx.strokeStyle = exp === -16 ? '#CBD5E1' : '#E2E8F0';
          stabCtx.beginPath();
          stabCtx.moveTo(padL, y); stabCtx.lineTo(leftW - 10, y);
          stabCtx.stroke();
          stabCtx.fillStyle = '#64748B';
          stabCtx.fillText(`1e${exp}`, padL - 5, y + 3);
        });

        const bars = [
          { label: 'CGS', val: eCGS.frob, col: '#DC2626' },
          { label: 'MGS', val: eMGS.frob, col: '#1D4ED8' },
          { label: 'Househ.', val: eHH.frob, col: '#0F766E' }
        ];
        const barSlot = (leftW - padL - 16) / 3;
        bars.forEach((b, idx) => {
          const logV = Math.max(-16.5, Math.min(0, Math.log10(Math.max(1e-17, b.val))));
          const fracFromBottom = (logV + 16.5) / 16.5;
          const bh = Math.max(4, fracFromBottom * chartH);
          const bx = padL + idx * barSlot + barSlot * 0.18;
          const bw = barSlot * 0.64;
          const by = padT + chartH - bh;

          stabCtx.fillStyle = b.col;
          stabCtx.fillRect(bx, by, bw, bh);

          stabCtx.fillStyle = '#0F172A';
          stabCtx.font = 'bold 10.5px sans-serif';
          stabCtx.textAlign = 'center';
          stabCtx.fillText(b.label, bx + bw / 2, h - 14);
          stabCtx.font = '9.5px monospace';
          stabCtx.fillStyle = b.col;
          stabCtx.fillText(fmtSci(b.val), bx + bw / 2, Math.max(padT + 10, by - 4));
        });

        // Right half: 3x3 Heatmap of |Q^T Q - I| for CGS vs Householder
        const hx0 = leftW + 12;
        const hw = w - hx0 - 10;
        stabCtx.fillStyle = '#0F172A';
        stabCtx.font = 'bold 11px sans-serif';
        stabCtx.textAlign = 'left';
        stabCtx.fillText('CGS |QᵀQ − I| Heatmap', hx0, 16);

        const cell = Math.min((hw - 6) / 3, (h - 54) / 3);
        const startY = 30;
        for (let r = 0; r < 3; r++) {
          for (let c = 0; c < 3; c++) {
            const err = eCGS.E[r][c];
            const intensity = Math.max(0, Math.min(1, (Math.log10(Math.max(1e-16, err)) + 16) / 16));
            const red = Math.round(248 - intensity * 30);
            const grn = Math.round(250 - intensity * 210);
            const blu = Math.round(252 - intensity * 210);
            const cx = hx0 + c * cell;
            const cy = startY + r * cell;
            stabCtx.fillStyle = `rgb(${red},${grn},${blu})`;
            stabCtx.fillRect(cx, cy, cell - 2, cell - 2);
            stabCtx.strokeStyle = '#CBD5E1';
            stabCtx.strokeRect(cx, cy, cell - 2, cell - 2);

            stabCtx.fillStyle = intensity > 0.65 ? '#FFFFFF' : '#0F172A';
            stabCtx.font = 'bold 9.5px monospace';
            stabCtx.textAlign = 'center';
            stabCtx.textBaseline = 'middle';
            const lbl = err < 1e-15 ? '0' : (err >= 0.01 ? err.toFixed(2) : `1e${Math.round(Math.log10(err))}`);
            stabCtx.fillText(lbl, cx + (cell - 2) / 2, cy + (cell - 2) / 2);
          }
        }
        stabCtx.textBaseline = 'alphabetic';
        stabCtx.fillStyle = '#64748B';
        stabCtx.font = '10px sans-serif';
        stabCtx.textAlign = 'left';
        stabCtx.fillText('Off-diagonals = |qᵢ · qⱼ| (ideal = 0)', hx0, startY + 3 * cell + 16);
      }

      epsSlider.addEventListener('input', () => {
        container.querySelectorAll('[data-eps]').forEach(b => b.classList.remove('active'));
        updateStabilityLab();
      });
      container.querySelectorAll('[data-eps]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-eps]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          epsSlider.value = btn.dataset.eps;
          updateStabilityLab();
        });
      });

      // =====================================================================
      // LAB 2: COSINE SIMILARITY & AI SEMANTIC SEARCH
      // =====================================================================
      const qSliders = [0, 1, 2, 3].map(i => container.querySelector(`#m9Q${i}`));
      const qVals = [0, 1, 2, 3].map(i => container.querySelector(`#m9Q${i}v`));
      const docScaleSlider = container.querySelector('#m9DocScale');
      const docScaleVal = container.querySelector('#m9DocScaleV');
      const cosCanvas = container.querySelector('#m9CosCanvas');
      const cosCtx = cosCanvas.getContext('2d');

      function updateCosineLab() {
        const q = qSliders.map((s, i) => {
          const v = parseFloat(s.value);
          qVals[i].textContent = v;
          return v;
        });
        const scaleB = parseFloat(docScaleSlider.value);
        docScaleVal.textContent = `${scaleB.toFixed(1)}×`;

        const nq = norm(q);
        const rows = DOCS.map(d => {
          const v = d.id === 'B' ? d.vec.map(x => x * scaleB) : d.vec.slice();
          const nv = norm(v);
          const cos = (nq > 0 && nv > 0) ? Math.max(-1, Math.min(1, dot(q, v) / (nq * nv))) : 0;
          const deg = (Math.acos(cos) * 180) / Math.PI;
          const dist = norm(q.map((x, i) => x - v[i]));
          return { ...d, scaledVec: v, normV: nv, cos, deg, dist };
        });

        // Compute rankings
        const byCos = rows.slice().sort((a, b) => b.cos - a.cos);
        const byEuc = rows.slice().sort((a, b) => a.dist - b.dist);
        rows.forEach(r => {
          r.cosRank = byCos.findIndex(x => x.id === r.id) + 1;
          r.eucRank = byEuc.findIndex(x => x.id === r.id) + 1;
        });

        // Sort display table by Cosine rank
        const sorted = rows.slice().sort((a, b) => a.cosRank - b.cosRank);
        const tbody = container.querySelector('#m9CosTableBody');
        tbody.innerHTML = sorted.map(r => {
          const highlight = r.id === 'B' ? 'background:#FFFBEB;' : '';
          const cosPill = r.cosRank === 1 ? 'top' : '';
          const eucPill = r.eucRank === 1 ? 'top' : (r.id === 'B' && r.eucRank > r.cosRank ? 'bad' : '');
          return `
            <tr style="${highlight}">
              <td><strong style="color:${r.color};">${r.name}</strong>${r.id === 'B' && scaleB !== 1 ? ` <span class="m9-pill warn">${scaleB.toFixed(1)}× len</span>` : ''}</td>
              <td style="font-family:monospace;">${r.deg.toFixed(1)}°</td>
              <td style="font-family:monospace;">${r.cos.toFixed(3)} ${cosPill ? `<span class="m9-pill ${cosPill}">#${r.cosRank}</span>` : `(#${r.cosRank})`}</td>
              <td style="font-family:monospace;">${r.dist.toFixed(2)} ${eucPill ? `<span class="m9-pill ${eucPill}">#${r.eucRank}</span>` : `(#${r.eucRank})`}</td>
            </tr>
          `;
        }).join('');

        const docB = rows.find(r => r.id === 'B');
        const insEl = container.querySelector('#m9CosInsight');
        if (byCos[0].id !== byEuc[0].id) {
          insEl.innerHTML = `
            <span class="m9-pill bad">Ranking Disagreement!</span>
            <strong>Cosine Search</strong> picks <strong>${byCos[0].name}</strong> ($\\theta = ${byCos[0].deg.toFixed(1)}^\\circ$), while <strong>Euclidean Distance</strong> picks <strong>${byEuc[0].name}</strong>!
            Notice Doc B is currently ranked <strong>#${docB.cosRank}</strong> by Cosine (invariant angle ${docB.deg.toFixed(1)}°) vs. <strong>#${docB.eucRank}</strong> by Euclidean distance ($\\|\\mathbf{q}-\\mathbf{d}_B\\| = ${docB.dist.toFixed(2)}$).
          `;
        } else {
          insEl.innerHTML = `
            <span class="m9-pill top">Rankings Agree (#1: ${byCos[0].id})</span>
            Try dragging the <strong>Doc B Length Scale ($s$)</strong> slider up to $2.5\\times$–$4.0\\times$: Doc B stretches outward along the <em>exact same ray</em>, keeping $\\cos\\theta = ${docB.cos.toFixed(3)}$ frozen while its Euclidean distance explodes!
          `;
        }
        renderMath(insEl);

        drawCosineCanvas(rows, nq);
      }

      function drawCosineCanvas(rows, nq) {
        if (!cosCanvas.isConnected) return;
        const dpr = window.devicePixelRatio || 1;
        const rect = cosCanvas.getBoundingClientRect();
        const w = rect.width || 380, h = rect.height || 220;
        cosCanvas.width = Math.round(w * dpr);
        cosCanvas.height = Math.round(h * dpr);
        cosCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
        cosCtx.clearRect(0, 0, w, h);

        const ox = 36, oy = h - 26;
        const maxR = Math.min(w - 64, h - 46);
        const unitR = maxR * 0.42; // represents norm = 14

        // Quarter-circle guide arcs (angle relative to query q along horizontal axis)
        [0.42, 0.84].forEach((f, idx) => {
          cosCtx.beginPath();
          cosCtx.arc(ox, oy, maxR * f, -Math.PI / 2, 0);
          cosCtx.strokeStyle = idx === 0 ? '#CBD5E1' : '#E2E8F0';
          cosCtx.setLineDash([4, 4]);
          cosCtx.stroke();
        });
        cosCtx.setLineDash([]);

        // Query ray along 0 degrees (horizontal)
        const qLen = Math.min(maxR, (nq / 14) * unitR);
        cosCtx.beginPath();
        cosCtx.moveTo(ox, oy);
        cosCtx.lineTo(ox + Math.max(24, qLen), oy);
        cosCtx.strokeStyle = '#0F172A';
        cosCtx.lineWidth = 3;
        cosCtx.stroke();
        cosCtx.beginPath();
        cosCtx.arc(ox + Math.max(24, qLen), oy, 5, 0, Math.PI * 2);
        cosCtx.fillStyle = '#0F172A';
        cosCtx.fill();
        cosCtx.font = 'bold 11px sans-serif';
        cosCtx.fillText('Query q (θ = 0°)', ox + Math.max(24, qLen) - 18, oy + 16);

        // Draw each document vector at its EXACT angle theta_i from q, with length proportional to ||d_i||
        rows.forEach(r => {
          const rad = (r.deg * Math.PI) / 180;
          const len = Math.min(maxR * 0.98, (r.normV / 14) * unitR);
          const ex = ox + len * Math.cos(rad);
          const ey = oy - len * Math.sin(rad);

          // Unit direction projection dot (on unit circle)
          const ux = ox + unitR * Math.cos(rad);
          const uy = oy - unitR * Math.sin(rad);

          cosCtx.beginPath();
          cosCtx.moveTo(ox, oy);
          cosCtx.lineTo(ex, ey);
          cosCtx.strokeStyle = r.color;
          cosCtx.lineWidth = r.id === 'B' ? 2.8 : 2;
          cosCtx.stroke();

          // Endpoint
          cosCtx.beginPath();
          cosCtx.arc(ex, ey, 4.5, 0, Math.PI * 2);
          cosCtx.fillStyle = r.color;
          cosCtx.fill();

          // Small hollow ring on normalized circle to show scale invariance
          cosCtx.beginPath();
          cosCtx.arc(ux, uy, 3, 0, Math.PI * 2);
          cosCtx.fillStyle = '#FFFFFF';
          cosCtx.fill();
          cosCtx.strokeStyle = r.color;
          cosCtx.lineWidth = 1.5;
          cosCtx.stroke();

          cosCtx.fillStyle = r.color;
          cosCtx.font = 'bold 10.5px sans-serif';
          cosCtx.fillText(`${r.id} (${r.deg.toFixed(0)}°)`, ex + 6, ey + 3);
        });

        cosCtx.fillStyle = '#475569';
        cosCtx.font = '10.5px sans-serif';
        cosCtx.fillText('Polar Angle Plot: angle = exact θ(q, dᵢ); length ∝ ||dᵢ||', 10, 15);
      }

      qSliders.forEach(s => s.addEventListener('input', () => {
        container.querySelectorAll('[data-qp]').forEach(b => b.classList.remove('active'));
        updateCosineLab();
      }));
      docScaleSlider.addEventListener('input', updateCosineLab);

      container.querySelectorAll('[data-qp]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-qp]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          const v = QUERY_PRESETS[btn.dataset.qp];
          v.forEach((val, i) => { qSliders[i].value = val; });
          updateCosineLab();
        });
      });

      // =====================================================================
      // LAB 3: ORTHONORMAL DCT COMPRESSION & SIGNAL DENOISING
      // =====================================================================
      let sigType = 'multitone';
      let noiseSeed = 2026;
      let unitNoise = makeGaussianNoise(noiseSeed);

      const kSlider = container.querySelector('#m9KSlider');
      const kValEl = container.querySelector('#m9KVal');
      const noiseSlider = container.querySelector('#m9NoiseSlider');
      const noiseValEl = container.querySelector('#m9NoiseVal');
      const sigCanvas = container.querySelector('#m9SigCanvas');
      const sigCtx = sigCanvas.getContext('2d');

      function updateSignalLab() {
        const K = parseInt(kSlider.value, 10);
        const sigma = parseFloat(noiseSlider.value);
        kValEl.textContent = `K = ${K} (${((K / SIG_N) * 100).toFixed(1)}%)`;
        noiseValEl.textContent = sigma.toFixed(2);

        const yClean = makeCleanSignal(sigType);
        const yNoisy = yClean.map((v, i) => v + sigma * unitNoise[i]);

        // Exact orthonormal DCT-II coordinates c_k = y_noisy . q_k
        const coeffs = DCT_BASIS.map(qk => dot(yNoisy, qk));

        // Orthogonal projection y_hat = sum_{k=0}^{K-1} c_k q_k
        const yHat = new Array(SIG_N).fill(0);
        for (let k = 0; k < K; k++) {
          const ck = coeffs[k];
          const qk = DCT_BASIS[k];
          for (let n = 0; n < SIG_N; n++) yHat[n] += ck * qk[n];
        }

        // Parseval & Denoising diagnostics
        const totalEnergy = dot(yNoisy, yNoisy);
        const parsevalSum = dot(coeffs, coeffs);
        let keptEnergy = 0;
        for (let k = 0; k < K; k++) keptEnergy += coeffs[k] * coeffs[k];
        const energyPct = totalEnergy > 0 ? (keptEnergy / totalEnergy) * 100 : 100;

        const noisyErr = norm(yNoisy.map((v, i) => v - yClean[i]));
        const denoisedErr = norm(yHat.map((v, i) => v - yClean[i]));
        const improved = denoisedErr < noisyErr;

        container.querySelector('#m9SigStats').innerHTML = `
          <div class="m9-stat-card">
            <div class="m9-stat-lbl">Parseval $\\|\\mathbf{y}\\|^2 = \\sum c_k^2$</div>
            <div class="m9-stat-val">${clean(totalEnergy, 2)} = ${clean(parsevalSum, 2)}</div>
          </div>
          <div class="m9-stat-card">
            <div class="m9-stat-lbl">Energy Retained ($W_K$)</div>
            <div class="m9-stat-val" style="color:#0F766E;">${energyPct.toFixed(1)}%</div>
          </div>
          <div class="m9-stat-card">
            <div class="m9-stat-lbl">Noisy Error $\\|\\mathbf{y}-\\mathbf{y}_{\\text{true}}\\|$</div>
            <div class="m9-stat-val" style="color:#64748B;">${clean(noisyErr, 3)}</div>
          </div>
          <div class="m9-stat-card">
            <div class="m9-stat-lbl">Projected Error $\\|\\hat{\\mathbf{y}}_K-\\mathbf{y}_{\\text{true}}\\|$</div>
            <div class="m9-stat-val" style="color:${improved ? '#166534' : '#DC2626'};">${clean(denoisedErr, 3)} ${improved ? '↓' : '↑'}</div>
          </div>
        `;
        renderMath(container.querySelector('#m9SigStats'));

        drawSignalCanvas(yClean, yNoisy, yHat, coeffs, K);
      }

      function drawSignalCanvas(yClean, yNoisy, yHat, coeffs, K) {
        if (!sigCanvas.isConnected) return;
        const dpr = window.devicePixelRatio || 1;
        const rect = sigCanvas.getBoundingClientRect();
        const w = rect.width || 700, h = rect.height || 290;
        sigCanvas.width = Math.round(w * dpr);
        sigCanvas.height = Math.round(h * dpr);
        sigCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
        sigCtx.clearRect(0, 0, w, h);

        const padL = 38, padR = 14;
        const plotW = w - padL - padR;

        // --- TOP PLOT: Time-domain signals (clean, noisy, projected) ---
        const topY = 20, topH = h * 0.52;
        const allY = yClean.concat(yNoisy, yHat);
        const yMin = Math.min(...allY) - 0.2;
        const yMax = Math.max(...allY) + 0.2;
        const mapX = i => padL + (i / (SIG_N - 1)) * plotW;
        const mapY = v => topY + topH - ((v - yMin) / (yMax - yMin || 1)) * topH;

        // Zero axis
        if (yMin < 0 && yMax > 0) {
          sigCtx.strokeStyle = '#E2E8F0';
          sigCtx.lineWidth = 1;
          sigCtx.beginPath();
          sigCtx.moveTo(padL, mapY(0)); sigCtx.lineTo(padL + plotW, mapY(0));
          sigCtx.stroke();
        }

        // Noisy signal (light slate)
        sigCtx.beginPath();
        yNoisy.forEach((v, i) => { if (i === 0) sigCtx.moveTo(mapX(i), mapY(v)); else sigCtx.lineTo(mapX(i), mapY(v)); });
        sigCtx.strokeStyle = 'rgba(148, 163, 184, 0.75)';
        sigCtx.lineWidth = 1.4;
        sigCtx.stroke();

        // Clean signal (dashed dark)
        sigCtx.beginPath();
        yClean.forEach((v, i) => { if (i === 0) sigCtx.moveTo(mapX(i), mapY(v)); else sigCtx.lineTo(mapX(i), mapY(v)); });
        sigCtx.strokeStyle = '#334155';
        sigCtx.lineWidth = 1.8;
        sigCtx.setLineDash([5, 4]);
        sigCtx.stroke();
        sigCtx.setLineDash([]);

        // Projected signal y_hat_K (bold teal)
        sigCtx.beginPath();
        yHat.forEach((v, i) => { if (i === 0) sigCtx.moveTo(mapX(i), mapY(v)); else sigCtx.lineTo(mapX(i), mapY(v)); });
        sigCtx.strokeStyle = '#0F766E';
        sigCtx.lineWidth = 2.6;
        sigCtx.stroke();

        // Top Legend
        sigCtx.font = 'bold 11px sans-serif';
        sigCtx.textAlign = 'left';
        sigCtx.fillStyle = '#0F766E'; sigCtx.fillText('━ Projection ŷ_K', padL, 13);
        sigCtx.fillStyle = '#334155'; sigCtx.fillText('┅ True Signal', padL + 115, 13);
        sigCtx.fillStyle = '#64748B'; sigCtx.fillText('━ Noisy Input y', padL + 215, 13);

        // --- BOTTOM PLOT: Energy Spectrum |c_k|^2 (Kept vs Discarded) ---
        const botY = topY + topH + 26, botH = h - botY - 18;
        const energies = coeffs.map(c => c * c);
        const maxE = Math.max(1e-6, ...energies);
        const bw = plotW / SIG_N;

        sigCtx.fillStyle = '#334155';
        sigCtx.font = 'bold 10.5px sans-serif';
        sigCtx.fillText(`DCT Energy Spectrum cₖ² (Teal = Kept modes 0..${K - 1}, Gray = Discarded noise modes ${K}..63)`, padL, botY - 6);

        for (let k = 0; k < SIG_N; k++) {
          // Mild sqrt compression so small noise bars remain visible
          const frac = Math.sqrt(energies[k] / maxE);
          const bh = Math.max(1.5, frac * botH);
          const bx = padL + k * bw;
          const by = botY + botH - bh;
          sigCtx.fillStyle = k < K ? '#0F766E' : '#CBD5E1';
          sigCtx.fillRect(bx + 0.6, by, Math.max(1, bw - 1.4), bh);
        }

        // Cutoff vertical line at K
        const cutX = padL + K * bw;
        sigCtx.strokeStyle = '#DC2626';
        sigCtx.lineWidth = 1.6;
        sigCtx.setLineDash([4, 3]);
        sigCtx.beginPath();
        sigCtx.moveTo(cutX, botY); sigCtx.lineTo(cutX, botY + botH);
        sigCtx.stroke();
        sigCtx.setLineDash([]);
      }

      kSlider.addEventListener('input', updateSignalLab);
      noiseSlider.addEventListener('input', updateSignalLab);
      container.querySelector('#m9ReseedBtn').addEventListener('click', () => {
        noiseSeed += 137;
        unitNoise = makeGaussianNoise(noiseSeed);
        updateSignalLab();
      });
      container.querySelectorAll('[data-sig]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-sig]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          sigType = btn.dataset.sig;
          updateSignalLab();
        });
      });

      // =====================================================================
      // SECTION D: CAPSTONE CONCEPTUAL CHECK
      // =====================================================================
      const quizWrap = container.querySelector('#m9CapstoneQuiz');
      quizWrap.innerHTML = CAPSTONE_QS.map((item, qi) => `
        <div class="m9-q-card" data-qcard="${qi}">
          <div style="font-size:0.89rem;"><strong>Q${qi + 1}.</strong> ${item.q}</div>
          ${item.opts.map((opt, oi) => `
            <button class="m9-q-opt" data-qi="${qi}" data-oi="${oi}">
              <strong>${String.fromCharCode(65 + oi)}.</strong> ${opt}
            </button>
          `).join('')}
          <div class="m9-q-exp" style="display:none; margin-top:8px; font-size:0.83rem;"></div>
        </div>
      `).join('');

      quizWrap.querySelectorAll('.m9-q-opt').forEach(btn => {
        btn.addEventListener('click', () => {
          const qi = parseInt(btn.dataset.qi, 10);
          const oi = parseInt(btn.dataset.oi, 10);
          const item = CAPSTONE_QS[qi];
          const card = quizWrap.querySelector(`[data-qcard="${qi}"]`);
          card.querySelectorAll('.m9-q-opt').forEach((b, idx) => {
            b.classList.remove('correct', 'wrong');
            if (idx === item.ans) b.classList.add('correct');
            else if (idx === oi) b.classList.add('wrong');
          });
          const expEl = card.querySelector('.m9-q-exp');
          expEl.style.display = 'block';
          expEl.style.color = oi === item.ans ? '#166534' : '#991B1B';
          expEl.innerHTML = `<strong>${oi === item.ans ? '✓ Correct!' : '✗ Not quite:'}</strong> ${item.exp}`;
          renderMath(expEl);
        });
      });
      renderMath(quizWrap);

      // Initial renders
      updateStabilityLab();
      updateCosineLab();
      updateSignalLab();

      window.addEventListener('resize', () => {
        updateStabilityLab();
        updateCosineLab();
        updateSignalLab();
      });
    }
  });
})();
