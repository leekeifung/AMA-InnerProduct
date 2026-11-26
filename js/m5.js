/**
 * AMA1751 Linear Algebra — Topic 6: Inner Product
 * Module 5: The Gram-Schmidt Process & QR Factorization
 */
(function () {
  // Self-contained vector & matrix utilities
  function dot(a, b) {
    let s = 0;
    for (let i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }
  function norm(a) { return Math.sqrt(dot(a, a)); }
  function scale(c, a) { return a.map(v => c * v); }
  function add(a, b) { return a.map((v, i) => v + b[i]); }
  function sub(a, b) { return a.map((v, i) => v - b[i]); }
  function transpose(M) {
    if (!M.length) return [];
    const r = M.length, c = M[0].length;
    return Array.from({ length: c }, (_, j) => Array.from({ length: r }, (_, i) => M[i][j]));
  }
  function matMul(A, B) {
    const r = A.length, k = A[0].length, c = B[0].length;
    const C = Array.from({ length: r }, () => Array(c).fill(0));
    for (let i = 0; i < r; i++) {
      for (let j = 0; j < c; j++) {
        let s = 0;
        for (let t = 0; t < k; t++) s += A[i][t] * B[t][j];
        C[i][j] = s;
      }
    }
    return C;
  }
  function maxAbsDiff(A, B) {
    let mx = 0;
    for (let i = 0; i < A.length; i++) {
      for (let j = 0; j < A[0].length; j++) {
        mx = Math.max(mx, Math.abs(A[i][j] - B[i][j]));
      }
    }
    return mx;
  }
  function eye(n) {
    return Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
  }
  function clean(x, d = 3) {
    if (!Number.isFinite(x)) return 'NaN';
    if (Math.abs(x) < 1e-11) return '0';
    if (Math.abs(x) < 1e-4) return x.toExponential(1);
    return String(Number(x.toFixed(d)));
  }
  function fmtVec(v, d = 2) {
    return `(${v.map(x => clean(x, d)).join(', ')})ᵀ`;
  }
  function renderMatHTML(M, d = 2) {
    let html = '<table class="m5-mat-table">';
    for (let i = 0; i < M.length; i++) {
      html += '<tr>';
      for (let j = 0; j < M[0].length; j++) {
        const val = M[i][j];
        const isZero = Math.abs(val) < 1e-11;
        html += `<td class="${isZero ? 'm5-zero-cell' : ''}">${clean(val, d)}</td>`;
      }
      html += '</tr>';
    }
    return html + '</table>';
  }

  /**
   * Detailed Gram-Schmidt on a list of column vectors `cols` (each length m).
   * mode: 'MGS' (Modified Gram-Schmidt) or 'CGS' (Classical Gram-Schmidt).
   */
  function runGramSchmidt(cols, mode = 'MGS') {
    const m = cols[0].length;
    const n = cols.length;
    const V = []; // orthogonal vectors (may contain zero vector if dependent)
    const Q = []; // unit vectors (or zero vector if dependent)
    const R = Array.from({ length: n }, () => Array(n).fill(0));
    const coeffsV = []; // classical projection ratios <x_j, v_i>/||v_i||^2
    const depFlags = [];

    for (let j = 0; j < n; j++) {
      const xj = cols[j].slice();
      let vj = xj.slice();
      const rowC = [];

      if (mode === 'CGS') {
        for (let i = 0; i < j; i++) {
          if (depFlags[i]) {
            rowC.push(0);
            continue;
          }
          const rij = dot(Q[i], xj);
          R[i][j] = rij;
          const vNormSq = dot(V[i], V[i]);
          rowC.push(vNormSq > 1e-24 ? dot(xj, V[i]) / vNormSq : 0);
          vj = sub(vj, scale(rij, Q[i]));
        }
      } else {
        // Modified Gram-Schmidt: update vj immediately against each prior q_i
        for (let i = 0; i < j; i++) {
          if (depFlags[i]) {
            rowC.push(0);
            continue;
          }
          const rij = dot(Q[i], vj);
          R[i][j] = rij;
          const vNormSq = dot(V[i], V[i]);
          rowC.push(vNormSq > 1e-24 ? dot(xj, V[i]) / vNormSq : 0);
          vj = sub(vj, scale(rij, Q[i]));
        }
      }

      const nv = norm(vj);
      const isDep = nv < 1e-10;
      depFlags.push(isDep);
      coeffsV.push(rowC);
      V.push(isDep ? Array(m).fill(0) : vj);
      R[j][j] = isDep ? 0 : nv;
      Q.push(isDep ? Array(m).fill(0) : scale(1 / nv, vj));
    }

    const Qmat = transpose(Q); // m x n
    const Amat = transpose(cols); // m x n
    const QtQ = matMul(transpose(Qmat), Qmat);
    const QR = matMul(Qmat, R);
    const hasDep = depFlags.some(Boolean);
    const orthoErr = hasDep ? NaN : maxAbsDiff(QtQ, eye(n));
    const reconErr = maxAbsDiff(QR, Amat);

    return { V, Q, Qmat, R, Amat, QtQ, QR, coeffsV, depFlags, hasDep, orthoErr, reconErr };
  }

  AMA1751.registerModule(5, {
    render(meta) {
      return `
        <style>
          .m5-canvas-wrap {
            position: relative;
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 10px;
            overflow: hidden;
          }
          .m5-canvas {
            display: block;
            width: 100%;
            height: 380px;
            cursor: grab;
            touch-action: none;
          }
          .m5-canvas:active { cursor: grabbing; }
          .m5-controls {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            margin-bottom: 12px;
            align-items: center;
          }
          .m5-chip {
            padding: 6px 11px;
            border-radius: 999px;
            border: 1px solid var(--border);
            background: #F8FAFC;
            font-size: 0.82rem;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.15s;
          }
          .m5-chip:hover, .m5-chip.active {
            background: var(--polyu-red-soft);
            border-color: #E8B4C4;
            color: var(--polyu-red-dark);
          }
          .m5-step-bar {
            display: flex;
            flex-wrap: wrap;
            gap: 6px;
            margin-bottom: 12px;
          }
          .m5-step-btn {
            flex: 1 1 120px;
            padding: 8px 10px;
            border-radius: 8px;
            border: 1px solid #CBD5E1;
            background: #FFFFFF;
            font-size: 0.82rem;
            font-weight: 700;
            cursor: pointer;
            color: #334155;
            transition: all 0.15s;
          }
          .m5-step-btn:hover { background: #F1F5F9; }
          .m5-step-btn.active {
            background: #0F766E;
            border-color: #0F766E;
            color: #FFFFFF;
            box-shadow: 0 2px 6px rgba(15, 118, 110, 0.25);
          }
          .m5-mat-table {
            border-collapse: separate;
            border-spacing: 3px;
            display: inline-table;
            padding: 5px 9px;
            border-left: 2px solid #475569;
            border-right: 2px solid #475569;
            border-radius: 6px;
            background: #FFFFFF;
          }
          .m5-mat-table td {
            min-width: 50px;
            text-align: center;
            padding: 4px 5px;
            border-radius: 4px;
            background: #F8FAFC;
            font-weight: 700;
            font-size: 0.83rem;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          }
          .m5-zero-cell { color: #94A3B8; font-weight: 500 !important; }
          .m5-plotly {
            width: 100%;
            height: 415px;
            border: 1px solid var(--border);
            border-radius: 10px;
            background: #FFFFFF;
            overflow: hidden;
          }
          .m5-log-box {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 12px 14px;
            font-size: 0.89rem;
            line-height: 1.62;
          }
          .m5-log-step {
            padding: 7px 9px;
            border-radius: 7px;
            margin-bottom: 7px;
            border-left: 3px solid #CBD5E1;
            background: #FFFFFF;
          }
          .m5-log-step.active-step {
            border-left-color: #0F766E;
            background: #F0FDFA;
          }
          .m5-qr-row {
            display: flex;
            flex-wrap: wrap;
            align-items: center;
            gap: 12px;
            padding: 12px;
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            overflow-x: auto;
          }
          .m5-qr-block {
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 4px;
          }
          .m5-qr-label {
            font-size: 0.76rem;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.04em;
            color: var(--text-muted);
          }
          .m5-eq-op {
            font-size: 1.25rem;
            font-weight: 800;
            color: #475569;
            padding-top: 16px;
          }
          .m5-edit-table {
            border-collapse: separate;
            border-spacing: 4px;
            padding: 6px 9px;
            border-left: 2px solid #1D4ED8;
            border-right: 2px solid #1D4ED8;
            border-radius: 6px;
            background: #FFFFFF;
          }
          .m5-edit-input {
            width: 64px;
            padding: 4px 5px;
            border: 1px solid #CBD5E1;
            border-radius: 5px;
            text-align: center;
            font-weight: 700;
            font-size: 0.85rem;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          }
          .m5-badge-row {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            margin-top: 10px;
          }
          .m5-badge {
            padding: 5px 10px;
            border-radius: 8px;
            font-size: 0.82rem;
            font-weight: 700;
            border: 1px solid transparent;
          }
          .m5-badge.ok { background: #DCFCE7; color: #166534; border-color: #86EFAC; }
          .m5-badge.warn { background: #FEF3C7; color: #92400E; border-color: #FCD34D; }
          .m5-badge.bad { background: #FEE2E2; color: #991B1B; border-color: #FCA5A5; }
          .m5-quiz-item {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 14px 16px;
            margin-bottom: 12px;
          }
          .m5-opt-group { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 10px; }
          .m5-opt-btn {
            padding: 7px 13px;
            border: 1px solid #CBD5E1;
            border-radius: 8px;
            background: #FFFFFF;
            cursor: pointer;
            font-weight: 600;
            font-size: 0.88rem;
          }
          .m5-opt-btn:hover { background: #F1F5F9; }
          .m5-opt-btn.correct { background: #DCFCE7; border-color: #16A34A; color: #166534; }
          .m5-opt-btn.wrong { background: #FEE2E2; border-color: #DC2626; color: #991B1B; }
          .m5-fb { margin-top: 8px; font-size: 0.88rem; font-weight: 500; }
        </style>

        <!-- HEADER -->
        <div class="module-header">
          <span class="slide-ref">Part 5 · Lecture Notes ${meta.slides}</span>
          <h2>5. The Gram-Schmidt Process &amp; $QR$ Factorization</h2>
          <p>How to manufacture an <strong>orthonormal basis</strong> from any linearly independent set by peeling off projections one vector at a time—and how the bookkeeping of that algorithm produces the matrix factorization $A = QR$.</p>
        </div>

        <!-- SECTION A: CORE ALGORITHM & THEOREMS -->
        <div class="card">
          <h3>The Gram-Schmidt Algorithm &amp; $QR$ Factorization Theorem</h3>
          <div class="grid-2">
            <div class="callout definition">
              <strong>Gram-Schmidt Process (from Basis $\\{\\mathbf{x}_1, \\dots, \\mathbf{x}_p\\}$ of $W$):</strong>
              $$\\begin{aligned}
                \\mathbf{v}_1 &= \\mathbf{x}_1 \\\\
                \\mathbf{v}_2 &= \\mathbf{x}_2 - \\frac{\\langle \\mathbf{x}_2, \\mathbf{v}_1 \\rangle}{\\|\\mathbf{v}_1\\|^2}\\mathbf{v}_1 \\\\
                \\mathbf{v}_k &= \\mathbf{x}_k - \\sum_{i=1}^{k-1} \\frac{\\langle \\mathbf{x}_k, \\mathbf{v}_i \\rangle}{\\|\\mathbf{v}_i\\|^2}\\mathbf{v}_i
              \\end{aligned}$$
              Finally, normalize: $\\mathbf{q}_k = \\dfrac{\\mathbf{v}_k}{\\|\\mathbf{v}_k\\|}$.<br/>
              <span style="font-size:0.86rem; color:var(--text-muted);"><strong>Nested Span Invariant:</strong> $\\operatorname{span}\\{\\mathbf{x}_1,\\dots,\\mathbf{x}_k\\} = \\operatorname{span}\\{\\mathbf{v}_1,\\dots,\\mathbf{v}_k\\} = \\operatorname{span}\\{\\mathbf{q}_1,\\dots,\\mathbf{q}_k\\}$ for every $k$.</span>
            </div>

            <div class="callout tip">
              <strong>Theorem ($QR$ Factorization):</strong><br/>
              If $A = [\\mathbf{x}_1 \\; \\cdots \\; \\mathbf{x}_n]$ is an $m \\times n$ matrix with <strong>linearly independent columns</strong>, then $A$ factors uniquely as
              $$A = QR$$
              <ul style="margin:6px 0 0 18px; padding:0; font-size:0.9rem; line-height:1.55;">
                <li>$Q = [\\mathbf{q}_1 \\; \\cdots \\; \\mathbf{q}_n]$ is $m \\times n$ with orthonormal columns ($Q^T Q = I_n$).</li>
                <li>$R = Q^T A$ is $n \\times n$ <strong>upper triangular</strong> with positive diagonal $r_{ii} = \\|\\mathbf{v}_i\\| > 0$.</li>
                <li><em>Why upper triangular?</em> Because $\\mathbf{x}_j \\in \\operatorname{span}\\{\\mathbf{q}_1,\\dots,\\mathbf{q}_j\\}$, so $r_{ij} = \\langle \\mathbf{q}_i, \\mathbf{x}_j \\rangle = 0$ for all $i > j$!</li>
              </ul>
            </div>
          </div>

          <div class="callout" style="margin-top:12px; background:#FFFBEB; border-left:4px solid #D97706;">
            <strong>Exam-Saving Tip (Avoid Messy Fractions!):</strong> You may rescale any intermediate orthogonal vector $\\mathbf{v}_k$ by a <em>positive</em> constant (e.g., multiply $(\\frac{1}{2}, -\\frac{1}{2}, 1)^T$ by $2$ to get $(1, -1, 2)^T$) before projecting the next vector—the direction of $\\mathbf{q}_k = \\mathbf{v}_k / \\|\\mathbf{v}_k\\|$ is unchanged! Then compute $R$ in one clean multiplication: <strong>$R = Q^T A$</strong>.
          </div>
        </div>

        <!-- SECTION B: INTERACTIVE LAB 1 — 2D GEOMETRIC STEPPER -->
        <div class="card">
          <h3>Interactive Lab 1: 2D Gram-Schmidt Step-by-Step Visualizer</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            Drag <strong style="color:#1D4ED8;">$\\mathbf{x}_1$ (blue)</strong> and <strong style="color:#7E22CE;">$\\mathbf{x}_2$ (purple)</strong> directly on the canvas, then step through the Gram-Schmidt stages to watch the shadow disappear and $A = QR$ assemble itself.
          </p>

          <div class="m5-controls">
            <strong style="font-size:0.84rem;">Presets:</strong>
            <button class="m5-chip active" data-p2="exam">Classic $(3,1)^T, (2,2)^T$</button>
            <button class="m5-chip" data-p2="skew">Very Skewed $(3,1)^T, (2.5,1.5)^T$</button>
            <button class="m5-chip" data-p2="obtuse">Obtuse Angle $(2,1)^T, (-2,2)^T$</button>
            <button class="m5-chip" data-p2="ortho">Already Orthogonal $(2,1)^T, (-1,2)^T$</button>
          </div>

          <div class="m5-step-bar" id="m5StepBar2D">
            <button class="m5-step-btn" data-step2="0">Step 0: Input $\\mathbf{x}_1, \\mathbf{x}_2$</button>
            <button class="m5-step-btn" data-step2="1">Step 1: Project $\\mathbf{x}_2$ onto $\\mathbf{v}_1$</button>
            <button class="m5-step-btn" data-step2="2">Step 2: Subtract $\\mathbf{v}_2 = \\mathbf{x}_2 - \\mathbf{p}$</button>
            <button class="m5-step-btn active" data-step2="3">Step 3: Normalize &amp; $A = QR$</button>
          </div>

          <div class="grid-2" style="align-items:start;">
            <div class="m5-canvas-wrap">
              <canvas id="m5Canvas2D" class="m5-canvas"></canvas>
            </div>

            <div>
              <div class="m5-log-box" id="m5Log2D"></div>

              <div class="m5-qr-row" style="margin-top:10px;" id="m5QR2D"></div>
            </div>
          </div>
        </div>

        <!-- SECTION C: INTERACTIVE LAB 2 — 3D GRAM-SCHMIDT STEPPER (PLOTLY 3D) -->
        <div class="card">
          <h3>Interactive Lab 2: 3D Gram-Schmidt Walkthrough &amp; Plane-Peeling Geometry</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            In Step 3, Gram-Schmidt subtracts the projection of $\\mathbf{x}_3$ onto the <strong>entire plane</strong> $W_2 = \\operatorname{span}\\{\\mathbf{v}_1, \\mathbf{v}_2\\}$ (shown in teal), leaving $\\mathbf{v}_3$ standing strictly perpendicular to $W_2$. Try the <em>Dependent Trap</em> preset to see what happens when $\\mathbf{x}_3 \\in W_2$!
          </p>

          <div class="m5-controls">
            <strong style="font-size:0.84rem;">3D Basis Preset:</strong>
            <button class="m5-chip active" data-p3="allInt">All-Integer Exam Basis</button>
            <button class="m5-chip" data-p3="sym">Symmetric $\\{(1,1,0)^T, (1,0,1)^T, (0,1,1)^T\\}$</button>
            <button class="m5-chip" data-p3="skew">Skewed Basis</button>
            <button class="m5-chip" data-p3="dep">Linearly Dependent Trap ($\\mathbf{v}_3 = \\mathbf{0}$)</button>
          </div>

          <div class="m5-step-bar" id="m5StepBar3D">
            <button class="m5-step-btn" data-step3="0">0. Input $\\mathbf{x}_1, \\mathbf{x}_2, \\mathbf{x}_3$</button>
            <button class="m5-step-btn" data-step3="1">1. Set $\\mathbf{v}_1 = \\mathbf{x}_1$</button>
            <button class="m5-step-btn" data-step3="2">2. Find $\\mathbf{v}_2 \\perp \\mathbf{v}_1$</button>
            <button class="m5-step-btn active" data-step3="3">3. Peel Plane $W_2 \\to \\mathbf{v}_3$</button>
            <button class="m5-step-btn" data-step3="4">4. Orthonormal Frame $\\mathbf{q}_1, \\mathbf{q}_2, \\mathbf{q}_3$</button>
          </div>

          <div class="grid-2" style="align-items:start;">
            <div id="m5Plotly3D" class="m5-plotly"></div>
            <div class="m5-log-box" id="m5Log3D"></div>
          </div>
        </div>

        <!-- SECTION D: LIVE A = QR MATRIX DECOMPOSER + CLASSICAL VS MODIFIED GS -->
        <div class="card">
          <h3>Interactive Lab 3: Live $A = QR$ Matrix Decomposer, Thin $QR$ &amp; Numerical Stability</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            Edit any entry of $A$ below (supports both square $3\\times 3$ and tall <strong>Thin $QR$</strong> $4\\times 3$ matrices) to compute $A = QR$ instantaneously. Try the <strong>Lauchli Ill-Conditioned</strong> preset to see why real numerical libraries use <em>Modified</em> Gram-Schmidt (or Householder) instead of <em>Classical</em> Gram-Schmidt!
          </p>

          <div class="m5-controls">
            <strong style="font-size:0.84rem;">Matrix Presets:</strong>
            <button class="m5-chip active" data-qrp="thin4x3">$4\\times 3$ Thin $QR$ (Exact Integer $R$!)</button>
            <button class="m5-chip" data-qrp="exam3x3">$3\\times 3$ Exam Matrix</button>
            <button class="m5-chip" data-qrp="dep3x3">$3\\times 3$ Rank-Deficient (Dependent Col 3)</button>
            <button class="m5-chip" data-qrp="lauchli">$4\\times 3$ Lauchli Ill-Conditioned ($\\varepsilon = 10^{-7}$)</button>
          </div>

          <div class="m5-controls" style="margin-bottom:10px;">
            <strong style="font-size:0.84rem;">Algorithm Mode:</strong>
            <button class="m5-chip active" data-gsmode="MGS">Modified Gram-Schmidt (MGS)</button>
            <button class="m5-chip" data-gsmode="CGS">Classical Gram-Schmidt (CGS)</button>
          </div>

          <div class="m5-qr-row">
            <div class="m5-qr-block">
              <div class="m5-qr-label" id="m5ALbl">Editable $A$</div>
              <div id="m5AEditor"></div>
            </div>
            <div class="m5-eq-op">=</div>
            <div class="m5-qr-block">
              <div class="m5-qr-label" id="m5QLbl">$Q$ (Orthonormal Cols)</div>
              <div id="m5QOut"></div>
            </div>
            <div class="m5-eq-op">·</div>
            <div class="m5-qr-block">
              <div class="m5-qr-label" id="m5RLbl">$R = Q^T A$ (Upper Tri.)</div>
              <div id="m5ROut"></div>
            </div>
          </div>

          <div class="m5-badge-row" id="m5QrBadges"></div>
          <div class="callout tip" id="m5QrInsight" style="margin-top:10px;"></div>

          <details class="code-details">
            <summary>Python (SymPy &amp; NumPy) Companion Code — Exact Gram-Schmidt &amp; $QR$ Factorization</summary>
            <pre><code>import sympy as sp
import numpy as np

# 1) Exact symbolic Gram-Schmidt & QR with SymPy (4x3 Thin QR)
A = sp.Matrix([
    [ 1, -1,  4],
    [ 1,  4, -2],
    [ 1,  4,  2],
    [ 1, -1,  0]
])

# Built-in exact QR decomposition (Q has orthonormal columns, R is upper triangular)
Q, R = A.QRdecomposition()
print("Exact Q =\\n", Q)
print("Exact R =\\n", R)
print("Verify Q^T Q == I:", Q.T * Q == sp.eye(3))   # True
print("Verify Q * R == A:", Q * R == A)             # True

# Step-by-step Gram-Schmidt on the columns of A
cols = [A.col(j) for j in range(A.cols)]
orthonormal_basis = sp.GramSchmidt(cols, orthonormal=True)

# 2) Fast numerical Thin QR with NumPy (uses Householder reflections internally)
A_np = np.array(A.tolist(), dtype=float)
Q_np, R_np = np.linalg.qr(A_np, mode='reduced')
# Note: NumPy may flip column signs of Q and rows of R; multiply by sign(diag(R))
# if you want strictly positive diagonal entries in R.
signs = np.sign(np.diag(R_np))
Q_pos = Q_np * signs
R_pos = signs[:, None] * R_np</code></pre>
          </details>
        </div>

        <!-- SECTION E: QUICK CONCEPT CHECK -->
        <div class="card">
          <h3>Quick Concept Check (Exam Pitfall Radar)</h3>

          <div class="m5-quiz-item" data-q="1" data-ans="B">
            <strong>Q1.</strong> In Step 3 of Gram-Schmidt, why must we compute
            $$\\mathbf{v}_3 = \\mathbf{x}_3 - \\frac{\\langle \\mathbf{x}_3, \\mathbf{v}_1\\rangle}{\\|\\mathbf{v}_1\\|^2}\\mathbf{v}_1 - \\frac{\\langle \\mathbf{x}_3, \\mathbf{v}_2\\rangle}{\\|\\mathbf{v}_2\\|^2}\\mathbf{v}_2$$
            using $\\mathbf{v}_1, \\mathbf{v}_2$ instead of the original vectors $\\mathbf{x}_1, \\mathbf{x}_2$?
            <div class="m5-opt-group">
              <button class="m5-opt-btn" data-opt="A">$\\mathbf{v}_1$ and $\\mathbf{v}_2$ always have integer entries, whereas $\\mathbf{x}_1, \\mathbf{x}_2$ do not.</button>
              <button class="m5-opt-btn" data-opt="B">The sum of 1D projections equals $\\operatorname{proj}_{W_2}(\\mathbf{x}_3)$ ONLY because $\\{\\mathbf{v}_1, \\mathbf{v}_2\\}$ is orthogonal.</button>
              <button class="m5-opt-btn" data-opt="C">Using $\\mathbf{x}_1, \\mathbf{x}_2$ would change the span of $W_2$.</button>
              <button class="m5-opt-btn" data-opt="D">Actually, using $\\mathbf{x}_1, \\mathbf{x}_2$ produces the exact same vector $\\mathbf{v}_3$.</button>
            </div>
            <div class="m5-fb" id="m5Fb1"></div>
          </div>

          <div class="m5-quiz-item" data-q="2" data-ans="C">
            <strong>Q2.</strong> While applying Gram-Schmidt to $\\{\\mathbf{x}_1, \\mathbf{x}_2, \\mathbf{x}_3\\}$, you obtain $\\mathbf{v}_1 \\ne \\mathbf{0}$, $\\mathbf{v}_2 \\ne \\mathbf{0}$, and $\\mathbf{v}_3 = \\mathbf{0}$. What does this tell you?
            <div class="m5-opt-group">
              <button class="m5-opt-btn" data-opt="A">You must have made an arithmetic error; Gram-Schmidt can never output $\\mathbf{0}$.</button>
              <button class="m5-opt-btn" data-opt="B">$\\mathbf{x}_3$ is already orthogonal to $\\mathbf{v}_1$ and $\\mathbf{v}_2$.</button>
              <button class="m5-opt-btn" data-opt="C">$\\mathbf{x}_3 \\in \\operatorname{span}\\{\\mathbf{x}_1, \\mathbf{x}_2\\}$ (linearly dependent!), so discard $\\mathbf{v}_3$—$\\{\\mathbf{v}_1, \\mathbf{v}_2\\}$ already spans the subspace.</button>
              <button class="m5-opt-btn" data-opt="D">The matrix $A = [\\mathbf{x}_1 \\; \\mathbf{x}_2 \\; \\mathbf{x}_3]$ is orthogonal.</button>
            </div>
            <div class="m5-fb" id="m5Fb2"></div>
          </div>

          <div class="m5-quiz-item" data-q="3" data-ans="A">
            <strong>Q3.</strong> Suppose $A$ has linearly independent columns and you have already constructed $Q$ via Gram-Schmidt so that $A = QR$ with $Q^T Q = I_n$. What is the fastest formula to compute $R$?
            <div class="m5-opt-group">
              <button class="m5-opt-btn" data-opt="A">$R = Q^T A$</button>
              <button class="m5-opt-btn" data-opt="B">$R = A Q^T$</button>
              <button class="m5-opt-btn" data-opt="C">$R = Q A$</button>
              <button class="m5-opt-btn" data-opt="D">$R = (A^T A)^{-1} Q$</button>
            </div>
            <div class="m5-fb" id="m5Fb3"></div>
          </div>

          <div class="m5-quiz-item" data-q="4" data-ans="D">
            <strong>Q4.</strong> Let $A$ be a $5 \\times 3$ matrix with linearly independent columns, and let $A = QR$ be its (thin) $QR$ factorization with positive diagonal entries in $R$. Which statement is <strong>TRUE</strong>?
            <div class="m5-opt-group">
              <button class="m5-opt-btn" data-opt="A">$Q$ is $5 \\times 5$ and $R$ is $5 \\times 3$, and $Q Q^T = I_5$.</button>
              <button class="m5-opt-btn" data-opt="B">$Q$ is $3 \\times 3$ and $R$ is $3 \\times 5$.</button>
              <button class="m5-opt-btn" data-opt="C">$Q$ is $5 \\times 3$, $R$ is $3 \\times 3$, and $Q Q^T = I_5$.</button>
              <button class="m5-opt-btn" data-opt="D">$Q$ is $5 \\times 3$ with $Q^T Q = I_3$, and $R$ is $3 \\times 3$ upper triangular with $r_{ii} = \\|\\mathbf{v}_i\\| > 0$.</button>
            </div>
            <div class="m5-fb" id="m5Fb4"></div>
          </div>
        </div>
      `;
    },

    mount(container, { renderMath }) {
      // =====================================================================
      // LAB 1: 2D GRAM-SCHMIDT GEOMETRIC STEPPER
      // =====================================================================
      const c2 = container.querySelector('#m5Canvas2D');
      const ctx2 = c2.getContext('2d');
      const presets2D = {
        exam:   { x1: [3, 1],   x2: [2, 2] },
        skew:   { x1: [3, 1],   x2: [2.5, 1.5] },
        obtuse: { x1: [2, 1],   x2: [-2, 2] },
        ortho:  { x1: [2, 1],   x2: [-1, 2] }
      };
      let x1_2d = presets2D.exam.x1.slice();
      let x2_2d = presets2D.exam.x2.slice();
      let step2D = 3;
      let drag2D = null; // 'x1' or 'x2'
      const range2D = 4.6;

      function layout2D() {
        const r = c2.getBoundingClientRect();
        const w = r.width || 450, h = r.height || 380;
        const sc = Math.min(w, h) / (2 * range2D);
        return { w, h, cx: w / 2, cy: h / 2, sc };
      }
      function toS2(p, L) { return [L.cx + p[0] * L.sc, L.cy - p[1] * L.sc]; }
      function toW2(sx, sy, L) {
        return [
          Math.max(-4.2, Math.min(4.2, (sx - L.cx) / L.sc)),
          Math.max(-4.2, Math.min(4.2, (L.cy - sy) / L.sc))
        ];
      }

      function drawArrow(ctx, x1, y1, x2, y2, color, width = 2.8, dashed = false) {
        const dx = x2 - x1, dy = y2 - y1;
        const len = Math.hypot(dx, dy);
        if (len < 2) return;
        const ang = Math.atan2(dy, dx), head = 9.5;
        ctx.save();
        ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = width;
        if (dashed) ctx.setLineDash([5, 4]);
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2 - head * 0.7 * Math.cos(ang), y2 - head * 0.7 * Math.sin(ang));
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.beginPath();
        ctx.moveTo(x2, y2);
        ctx.lineTo(x2 - head * Math.cos(ang - Math.PI / 7), y2 - head * Math.sin(ang - Math.PI / 7));
        ctx.lineTo(x2 - head * Math.cos(ang + Math.PI / 7), y2 - head * Math.sin(ang + Math.PI / 7));
        ctx.closePath();
        ctx.fill();
        ctx.restore();
      }

      function drawLab2D() {
        if (!c2.isConnected) return;
        const dpr = window.devicePixelRatio || 1;
        const L = layout2D();
        c2.width = Math.round(L.w * dpr);
        c2.height = Math.round(L.h * dpr);
        ctx2.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx2.clearRect(0, 0, L.w, L.h);

        // Grid & axes
        ctx2.strokeStyle = '#F1F5F9'; ctx2.lineWidth = 1;
        for (let g = -4; g <= 4; g++) {
          const [gx] = toS2([g, 0], L), [, gy] = toS2([0, g], L);
          ctx2.beginPath(); ctx2.moveTo(gx, 0); ctx2.lineTo(gx, L.h); ctx2.stroke();
          ctx2.beginPath(); ctx2.moveTo(0, gy); ctx2.lineTo(L.w, gy); ctx2.stroke();
        }
        ctx2.strokeStyle = '#CBD5E1'; ctx2.lineWidth = 1.2;
        ctx2.beginPath(); ctx2.moveTo(0, L.cy); ctx2.lineTo(L.w, L.cy); ctx2.stroke();
        ctx2.beginPath(); ctx2.moveTo(L.cx, 0); ctx2.lineTo(L.cx, L.h); ctx2.stroke();

        const n1 = norm(x1_2d);
        const gs = runGramSchmidt([x1_2d, x2_2d], 'MGS');
        const v1 = gs.V[0], v2 = gs.V[1];
        const q1 = gs.Q[0], q2 = gs.Q[1];
        const mu = n1 > 1e-8 ? dot(x2_2d, v1) / (n1 * n1) : 0;
        const proj = scale(mu, v1);

        // Span(v1) faint guide line in Steps >= 1
        if (step2D >= 1 && n1 > 1e-6) {
          const u1 = scale(1 / n1, v1);
          const a = toS2(scale(-7, u1), L), b = toS2(scale(7, u1), L);
          ctx2.save();
          ctx2.strokeStyle = 'rgba(29, 78, 216, 0.25)'; ctx2.lineWidth = 2; ctx2.setLineDash([6, 4]);
          ctx2.beginPath(); ctx2.moveTo(a[0], a[1]); ctx2.lineTo(b[0], b[1]); ctx2.stroke();
          ctx2.restore();
        }

        // Unit circle in Step 3
        if (step2D === 3) {
          ctx2.save();
          ctx2.strokeStyle = 'rgba(15, 118, 110, 0.35)';
          ctx2.fillStyle = 'rgba(15, 118, 110, 0.04)';
          ctx2.lineWidth = 1.6; ctx2.setLineDash([5, 4]);
          ctx2.beginPath();
          ctx2.arc(L.cx, L.cy, L.sc, 0, Math.PI * 2);
          ctx2.fill(); ctx2.stroke();
          ctx2.restore();
        }

        const [x1x, x1y] = toS2(x1_2d, L);
        const [x2x, x2y] = toS2(x2_2d, L);
        const [px, py] = toS2(proj, L);
        const [v2x, v2y] = toS2(v2, L);

        // Draw original x1, x2 (slightly dimmed in Step 3 so q1, q2 pop)
        const dimAlpha = step2D === 3 ? 0.38 : 1.0;
        ctx2.save();
        ctx2.globalAlpha = dimAlpha;
        drawArrow(ctx2, L.cx, L.cy, x1x, x1y, '#1D4ED8', 3.0);
        drawArrow(ctx2, L.cx, L.cy, x2x, x2y, '#7E22CE', 3.0);
        ctx2.restore();

        // Step 1+: Projection shadow p and dashed perpendicular drop from p to x2
        if (step2D >= 1 && n1 > 1e-6) {
          ctx2.save();
          ctx2.globalAlpha = step2D === 3 ? 0.35 : 0.95;
          drawArrow(ctx2, L.cx, L.cy, px, py, '#D97706', 2.8);
          ctx2.strokeStyle = '#0F766E'; ctx2.lineWidth = 2.0; ctx2.setLineDash([5, 4]);
          ctx2.beginPath(); ctx2.moveTo(px, py); ctx2.lineTo(x2x, x2y); ctx2.stroke();
          ctx2.restore();

          if (step2D <= 2) {
            ctx2.fillStyle = '#B45309'; ctx2.font = 'bold 12px sans-serif';
            ctx2.fillText('p = proj(x₂)', px + 7, py + 15);
          }
        }

        // Step 2+: Orthogonal vector v2 at origin + right-angle square
        if (step2D >= 2 && !gs.hasDep) {
          ctx2.save();
          ctx2.globalAlpha = step2D === 3 ? 0.35 : 1.0;
          drawArrow(ctx2, L.cx, L.cy, v2x, v2y, '#0F766E', 3.2);
          ctx2.restore();

          // Right-angle marker at origin between q1 and q2
          const sq = 0.28;
          const p1 = toS2(scale(sq, q1), L);
          const p2 = toS2(add(scale(sq, q1), scale(sq, q2)), L);
          const p3 = toS2(scale(sq, q2), L);
          ctx2.save();
          ctx2.strokeStyle = '#0F766E'; ctx2.lineWidth = 1.6;
          ctx2.beginPath();
          ctx2.moveTo(p1[0], p1[1]); ctx2.lineTo(p2[0], p2[1]); ctx2.lineTo(p3[0], p3[1]);
          ctx2.stroke();
          ctx2.restore();

          if (step2D === 2) {
            ctx2.fillStyle = '#0F766E'; ctx2.font = 'bold 13px sans-serif';
            ctx2.fillText('v₂ ⊥ v₁', v2x + 8, v2y - 6);
          }
        }

        // Step 3: Orthonormal unit vectors q1, q2
        if (step2D === 3 && !gs.hasDep) {
          const [q1x, q1y] = toS2(q1, L);
          const [q2x, q2y] = toS2(q2, L);
          drawArrow(ctx2, L.cx, L.cy, q1x, q1y, '#B91C1C', 3.4);
          drawArrow(ctx2, L.cx, L.cy, q2x, q2y, '#0F766E', 3.4);
          ctx2.font = 'bold 14px sans-serif';
          ctx2.fillStyle = '#B91C1C'; ctx2.fillText('q₁', q1x + 8, q1y - 6);
          ctx2.fillStyle = '#0F766E'; ctx2.fillText('q₂', q2x + 8, q2y - 6);
        }

        // Draggable handles for x1 and x2
        [
          { sx: x1x, sy: x1y, col: '#1D4ED8', lbl: 'x₁' },
          { sx: x2x, sy: x2y, col: '#7E22CE', lbl: 'x₂' }
        ].forEach(h => {
          ctx2.beginPath(); ctx2.arc(h.sx, h.sy, 6.8, 0, Math.PI * 2);
          ctx2.fillStyle = '#FFFFFF'; ctx2.fill();
          ctx2.lineWidth = 2.4; ctx2.strokeStyle = h.col; ctx2.stroke();
          ctx2.fillStyle = h.col; ctx2.font = 'bold 13px sans-serif';
          ctx2.fillText(h.lbl, h.sx + 9, h.sy - 7);
        });

        // Update algebraic log
        const logEl = container.querySelector('#m5Log2D');
        const qrEl = container.querySelector('#m5QR2D');
        if (gs.hasDep) {
          logEl.innerHTML = `
            <div class="m5-log-step active-step" style="border-left-color:#DC2626; background:#FEF2F2;">
              <strong style="color:#991B1B;">Linearly Dependent Vectors Detected!</strong><br/>
              $\\mathbf{x}_1 = ${fmtVec(x1_2d)}$ and $\\mathbf{x}_2 = ${fmtVec(x2_2d)}$ are collinear, so $\\mathbf{v}_2 = \\mathbf{0}$. Drag either vector off the line to form a basis of $\\mathbb{R}^2$.
            </div>`;
          qrEl.innerHTML = `<span style="color:#991B1B; font-weight:700; font-size:0.86rem;">$A = QR$ with invertible $R$ requires linearly independent columns.</span>`;
          renderMath(logEl); renderMath(qrEl);
          return;
        }

        logEl.innerHTML = `
          <div class="m5-log-step ${step2D === 0 ? 'active-step' : ''}">
            <strong>Step 0 (Input):</strong> $\\mathbf{x}_1 = ${fmtVec(x1_2d)}$, $\\mathbf{x}_2 = ${fmtVec(x2_2d)}$
          </div>
          <div class="m5-log-step ${step2D === 1 ? 'active-step' : ''}">
            <strong>Step 1 (Project onto $\\mathbf{v}_1 = \\mathbf{x}_1$):</strong>
            $\\mathbf{p} = \\dfrac{\\langle \\mathbf{x}_2, \\mathbf{v}_1 \\rangle}{\\|\\mathbf{v}_1\\|^2}\\mathbf{v}_1 = \\dfrac{${clean(dot(x2_2d, v1), 2)}}{${clean(n1 * n1, 2)}}\\mathbf{v}_1 = ${fmtVec(proj)}$
          </div>
          <div class="m5-log-step ${step2D === 2 ? 'active-step' : ''}">
            <strong>Step 2 (Peel Shadow):</strong>
            $\\mathbf{v}_2 = \\mathbf{x}_2 - \\mathbf{p} = ${fmtVec(v2)}$ &nbsp;(check $\\langle \\mathbf{v}_1, \\mathbf{v}_2\\rangle = ${clean(dot(v1, v2), 3)}$ ✓)
          </div>
          <div class="m5-log-step ${step2D === 3 ? 'active-step' : ''}">
            <strong>Step 3 (Normalize):</strong>
            $\\mathbf{q}_1 = \\dfrac{\\mathbf{v}_1}{${clean(norm(v1), 2)}} = ${fmtVec(q1)}$, &nbsp;
            $\\mathbf{q}_2 = \\dfrac{\\mathbf{v}_2}{${clean(norm(v2), 2)}} = ${fmtVec(q2)}$
          </div>
        `;
        renderMath(logEl);

        qrEl.innerHTML = `
          <div class="m5-qr-block"><div class="m5-qr-label">$A$</div>${renderMatHTML(gs.Amat, 2)}</div>
          <div class="m5-eq-op">=</div>
          <div class="m5-qr-block"><div class="m5-qr-label">$Q$</div>${renderMatHTML(gs.Qmat, 2)}</div>
          <div class="m5-eq-op">·</div>
          <div class="m5-qr-block"><div class="m5-qr-label">$R = Q^T A$</div>${renderMatHTML(gs.R, 2)}</div>
        `;
        renderMath(qrEl);
      }

      function ptr2(e) {
        const r = c2.getBoundingClientRect();
        const cx = e.touches ? e.touches[0].clientX : e.clientX;
        const cy = e.touches ? e.touches[0].clientY : e.clientY;
        return [cx - r.left, cy - r.top];
      }
      function startDrag2(e) {
        const L = layout2D(), [px, py] = ptr2(e);
        const [x1x, x1y] = toS2(x1_2d, L), [x2x, x2y] = toS2(x2_2d, L);
        if (Math.hypot(px - x1x, py - x1y) < 24) { drag2D = 'x1'; e.preventDefault(); }
        else if (Math.hypot(px - x2x, py - x2y) < 24) { drag2D = 'x2'; e.preventDefault(); }
      }
      function moveDrag2(e) {
        if (!drag2D || !c2.isConnected) return;
        if (e.touches) e.preventDefault();
        const L = layout2D(), [px, py] = ptr2(e), pt = toW2(px, py, L);
        const snapped = [Math.round(pt[0] * 4) / 4, Math.round(pt[1] * 4) / 4];
        if (drag2D === 'x1') x1_2d = snapped;
        else x2_2d = snapped;
        drawLab2D();
      }
      c2.addEventListener('mousedown', startDrag2);
      c2.addEventListener('touchstart', startDrag2, { passive: false });
      window.addEventListener('mousemove', moveDrag2);
      window.addEventListener('touchmove', moveDrag2, { passive: false });
      window.addEventListener('mouseup', () => { drag2D = null; });
      window.addEventListener('touchend', () => { drag2D = null; });

      container.querySelectorAll('[data-p2]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-p2]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          const p = presets2D[btn.dataset.p2];
          x1_2d = p.x1.slice(); x2_2d = p.x2.slice();
          drawLab2D();
        });
      });
      container.querySelectorAll('[data-step2]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-step2]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          step2D = parseInt(btn.dataset.step2, 10);
          drawLab2D();
        });
      });

      // =====================================================================
      // LAB 2: 3D GRAM-SCHMIDT WALKTHROUGH (PLOTLY 3D)
      // =====================================================================
      const presets3D = {
        allInt: { cols: [[1, 1, 1], [1, 0, 2], [0, 2, 4]], label: 'All-Integer Orthogonal Steps' },
        sym:    { cols: [[1, 1, 0], [1, 0, 1], [0, 1, 1]], label: 'Symmetric Basis' },
        skew:   { cols: [[2, 1, 0], [1, 2, 1], [0, 1, 2]], label: 'General Skewed Basis' },
        dep:    { cols: [[1, 1, 0], [0, 1, 1], [2, 3, 1]], label: 'Linearly Dependent Trap (x₃ = 2x₁ + x₂)' }
      };
      let activeP3 = 'allInt';
      let step3D = 3;
      const p3El = container.querySelector('#m5Plotly3D');

      function vecTrace(v, name, color, label, width = 6, dash = 'solid', opacity = 1) {
        return {
          type: 'scatter3d', mode: 'lines+markers+text', name,
          x: [0, v[0]], y: [0, v[1]], z: [0, v[2]],
          text: ['', label], textposition: 'top center',
          textfont: { color, size: 13 },
          opacity,
          line: { color, width, dash },
          marker: { size: [2, 5.5], color }
        };
      }

      function updateLab3D() {
        const cols = presets3D[activeP3].cols;
        const [x1, x2, x3] = cols;
        const gs = runGramSchmidt(cols, 'MGS');
        const [v1, v2, v3] = gs.V;
        const [q1, q2, q3] = gs.Q;
        const c21 = gs.coeffsV[1][0];
        const c31 = gs.coeffsV[2][0];
        const c32 = gs.coeffsV[2][1];
        const p2 = scale(c21, v1);
        const p3 = add(scale(c31, v1), scale(c32, v2));

        // Build algebraic step log
        const log3 = container.querySelector('#m5Log3D');
        const depNote = gs.depFlags[2]
          ? `<br/><span style="color:#DC2626; font-weight:700;">⚠ $\\mathbf{v}_3 = \\mathbf{0}$! Notice $\\mathbf{x}_3 = \\operatorname{proj}_{W_2}(\\mathbf{x}_3)$, proving $\\mathbf{x}_3$ lies completely inside $\\operatorname{span}\\{\\mathbf{v}_1, \\mathbf{v}_2\\}$ (linearly dependent!). Discard $\\mathbf{v}_3$.</span>`
          : `<br/><span style="color:#15803D;">Check: $\\langle \\mathbf{v}_3, \\mathbf{v}_1\\rangle = ${clean(dot(v3, v1), 3)}$, $\\langle \\mathbf{v}_3, \\mathbf{v}_2\\rangle = ${clean(dot(v3, v2), 3)}$ ✓</span>`;

        log3.innerHTML = `
          <div class="m5-log-step ${step3D === 0 ? 'active-step' : ''}">
            <strong>Step 0 (Input Columns):</strong><br/>
            $\\mathbf{x}_1 = ${fmtVec(x1)}$, $\\mathbf{x}_2 = ${fmtVec(x2)}$, $\\mathbf{x}_3 = ${fmtVec(x3)}$
          </div>
          <div class="m5-log-step ${step3D === 1 ? 'active-step' : ''}">
            <strong>Step 1 (Anchor First Direction):</strong><br/>
            $\\mathbf{v}_1 = \\mathbf{x}_1 = ${fmtVec(v1)}$, &nbsp;$\\|\\mathbf{v}_1\\|^2 = ${clean(dot(v1, v1), 2)}$
          </div>
          <div class="m5-log-step ${step3D === 2 ? 'active-step' : ''}">
            <strong>Step 2 (Orthogonalize $\\mathbf{x}_2$ against $\\mathbf{v}_1$):</strong><br/>
            $\\mathbf{v}_2 = \\mathbf{x}_2 - (${clean(c21, 2)})\\mathbf{v}_1 = ${fmtVec(v2)}$, &nbsp;$\\|\\mathbf{v}_2\\|^2 = ${clean(dot(v2, v2), 2)}$
          </div>
          <div class="m5-log-step ${step3D === 3 ? 'active-step' : ''}">
            <strong>Step 3 (Peel Plane $W_2 = \\operatorname{span}\\{\\mathbf{v}_1, \\mathbf{v}_2\\}$):</strong><br/>
            $\\operatorname{proj}_{W_2}(\\mathbf{x}_3) = (${clean(c31, 2)})\\mathbf{v}_1 + (${clean(c32, 2)})\\mathbf{v}_2 = ${fmtVec(p3)}$<br/>
            $\\mathbf{v}_3 = \\mathbf{x}_3 - \\operatorname{proj}_{W_2}(\\mathbf{x}_3) = ${fmtVec(v3)}$${depNote}
          </div>
          <div class="m5-log-step ${step3D === 4 ? 'active-step' : ''}">
            <strong>Step 4 (Normalize to Orthonormal Frame):</strong><br/>
            $\\mathbf{q}_1 = ${fmtVec(q1)}$, $\\mathbf{q}_2 = ${fmtVec(q2)}$${gs.depFlags[2] ? ', (2D subspace basis complete!)' : `, $\\mathbf{q}_3 = ${fmtVec(q3)}$`}
          </div>
        `;
        renderMath(log3);

        if (typeof Plotly === 'undefined') {
          p3El.innerHTML = '<div style="padding:20px; color:var(--text-muted);">3D plot requires Plotly.js.</div>';
          return;
        }

        const traces = [];
        // Input vectors x1, x2, x3
        const origOpacity = step3D === 0 ? 1.0 : 0.35;
        traces.push(vecTrace(x1, 'x₁', '#1D4ED8', 'x₁', step3D === 0 ? 6 : 3, step3D === 0 ? 'solid' : 'dot', origOpacity));
        traces.push(vecTrace(x2, 'x₂', '#7E22CE', 'x₂', step3D === 0 ? 6 : 3, step3D === 0 ? 'solid' : 'dot', origOpacity));
        traces.push(vecTrace(x3, 'x₃', '#DC2626', 'x₃', step3D === 0 || step3D === 3 ? 5.5 : 3, step3D === 0 || step3D === 3 ? 'solid' : 'dot', step3D === 0 || step3D === 3 ? 0.95 : 0.3));

        // Step 1+: v1
        if (step3D >= 1 && step3D <= 3) {
          traces.push(vecTrace(v1, 'v₁', '#1D4ED8', 'v₁', 6.5));
        }
        // Step 2: projection of x2 onto v1 and v2
        if (step3D === 2) {
          traces.push(vecTrace(p2, 'proj_{v₁}(x₂)', '#D97706', 'p₂', 5, 'dash'));
          traces.push({
            type: 'scatter3d', mode: 'lines', showlegend: false,
            x: [p2[0], x2[0]], y: [p2[1], x2[1]], z: [p2[2], x2[2]],
            line: { color: '#0F766E', width: 4, dash: 'dash' }
          });
        }
        if (step3D >= 2 && step3D <= 3 && !gs.depFlags[1]) {
          traces.push(vecTrace(v2, 'v₂', '#0F766E', 'v₂', 6.5));
        }

        // Step 3: Plane W2 = span{v1, v2}, shadow p3, drop segment, and v3
        if (step3D === 3 && !gs.depFlags[1]) {
          const gv = [-3.4, 0, 3.4];
          const sx = [], sy = [], sz = [];
          for (const u of gv) {
            const rx = [], ry = [], rz = [];
            for (const w of gv) {
              const pt = add(scale(u, q1), scale(w, q2));
              rx.push(pt[0]); ry.push(pt[1]); rz.push(pt[2]);
            }
            sx.push(rx); sy.push(ry); sz.push(rz);
          }
          traces.push({
            type: 'surface', name: 'Plane W₂ = span{v₁,v₂}', showscale: false, opacity: 0.28,
            x: sx, y: sy, z: sz,
            colorscale: [[0, '#99F6E4'], [1, '#0F766E']], hoverinfo: 'skip'
          });
          traces.push(vecTrace(p3, 'proj_{W₂}(x₃)', '#D97706', 'p₃', 5, 'dash'));
          if (!gs.depFlags[2]) {
            traces.push({
              type: 'scatter3d', mode: 'lines', showlegend: false,
              x: [p3[0], x3[0]], y: [p3[1], x3[1]], z: [p3[2], x3[2]],
              line: { color: '#B91C1C', width: 4.5, dash: 'dash' }
            });
            traces.push(vecTrace(v3, 'v₃ ⊥ W₂', '#B91C1C', 'v₃', 7));
          }
        }

        // Step 4: Orthonormal frame q1, q2, q3
        if (step3D === 4) {
          traces.push(vecTrace(q1, 'q₁ (unit)', '#1D4ED8', 'q₁', 7));
          if (!gs.depFlags[1]) traces.push(vecTrace(q2, 'q₂ (unit)', '#0F766E', 'q₂', 7));
          if (!gs.depFlags[2]) traces.push(vecTrace(q3, 'q₃ (unit)', '#B91C1C', 'q₃', 7));
        }

        const lim = activeP3 === 'allInt' ? 4.4 : 3.4;
        Plotly.react(p3El, traces, {
          margin: { l: 0, r: 0, t: 8, b: 0 },
          legend: { orientation: 'h', y: 0.99, x: 0.01, font: { size: 10.5 } },
          scene: {
            aspectmode: 'cube',
            xaxis: { range: [-lim, lim], title: 'x₁', zerolinecolor: '#94A3B8' },
            yaxis: { range: [-lim, lim], title: 'x₂', zerolinecolor: '#94A3B8' },
            zaxis: { range: [-lim, lim], title: 'x₃', zerolinecolor: '#94A3B8' },
            camera: { eye: { x: 1.55, y: 1.45, z: 1.15 } }
          }
        }, { responsive: true, displayModeBar: false });
      }

      container.querySelectorAll('[data-p3]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-p3]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          activeP3 = btn.dataset.p3;
          updateLab3D();
        });
      });
      container.querySelectorAll('[data-step3]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-step3]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          step3D = parseInt(btn.dataset.step3, 10);
          updateLab3D();
        });
      });

      // =====================================================================
      // LAB 3: LIVE A = QR MATRIX DECOMPOSER + CGS VS MGS
      // =====================================================================
      const qrPresets = {
        thin4x3: {
          rows: 4, cols: 3,
          mat: [
            [1, -1, 4],
            [1,  4, -2],
            [1,  4, 2],
            [1, -1, 0]
          ],
          note: '<strong>Thin $QR$ Gem ($4\\times 3$):</strong> Every entry of $Q$ is $\\pm 0.5$ and every entry of $R$ is an exact integer! Notice $Q^T Q = I_3$, while $Q Q^T = P_{\\operatorname{Col}(A)} \\ne I_4$ is the $4\\times 4$ orthogonal projection matrix onto $\\operatorname{Col}(A)$.'
        },
        exam3x3: {
          rows: 3, cols: 3,
          mat: [
            [1, 1, 0],
            [1, 0, 2],
            [1, 2, 4]
          ],
          note: '<strong>All-Integer Orthogonal Steps:</strong> Diagonal entries of $R$ are $r_{11} = \\sqrt{3} \\approx 1.73$, $r_{22} = \\sqrt{2} \\approx 1.41$, and $r_{33} = \\sqrt{6} \\approx 2.45$, equal to $\\|\\mathbf{v}_1\\|, \\|\\mathbf{v}_2\\|, \\|\\mathbf{v}_3\\|$.'
        },
        dep3x3: {
          rows: 3, cols: 3,
          mat: [
            [1, 0, 2],
            [1, 1, 3],
            [0, 1, 1]
          ],
          note: '<strong>Rank-Deficient Warning:</strong> Column 3 equals $2\\mathbf{x}_1 + \\mathbf{x}_2$, so $r_{33} = \\|\\mathbf{v}_3\\| = 0$! Full $QR$ with an invertible upper-triangular $R$ requires $A$ to have linearly independent columns.'
        },
        lauchli: {
          rows: 4, cols: 3,
          mat: [
            [1, 1, 1],
            [1e-7, 0, 0],
            [0, 1e-7, 0],
            [0, 0, 1e-7]
          ],
          note: '<strong>Numerical Stability Classic (Lauchli Matrix, $\\varepsilon = 10^{-7}$):</strong> Toggle between <strong>MGS</strong> and <strong>CGS</strong> above! In floating-point arithmetic, Classical GS loses orthogonality catastrophically ($\\|Q^TQ - I\\|_{\\max} \\approx 0.50$), whereas Modified GS maintains $\\approx 10^{-9}$ by re-projecting the updated residual at each step.'
        }
      };

      let activeQrKey = 'thin4x3';
      let gsMode = 'MGS';
      let currentMat = qrPresets.thin4x3.mat.map(r => r.slice());

      function buildAEditor() {
        const m = currentMat.length, n = currentMat[0].length;
        let html = '<table class="m5-edit-table">';
        for (let i = 0; i < m; i++) {
          html += '<tr>';
          for (let j = 0; j < n; j++) {
            const valStr = Math.abs(currentMat[i][j]) > 0 && Math.abs(currentMat[i][j]) < 1e-4
              ? currentMat[i][j].toExponential(0)
              : String(currentMat[i][j]);
            html += `<td><input type="number" step="any" class="m5-edit-input" data-r="${i}" data-c="${j}" value="${valStr}" /></td>`;
          }
          html += '</tr>';
        }
        html += '</table>';
        const ed = container.querySelector('#m5AEditor');
        ed.innerHTML = html;
        ed.querySelectorAll('.m5-edit-input').forEach(inp => {
          inp.addEventListener('input', () => {
            const r = parseInt(inp.dataset.r, 10);
            const c = parseInt(inp.dataset.c, 10);
            const v = parseFloat(inp.value);
            currentMat[r][c] = Number.isFinite(v) ? v : 0;
            computeQRDisplay();
          });
        });
      }

      function computeQRDisplay() {
        const m = currentMat.length, n = currentMat[0].length;
        const cols = transpose(currentMat);
        const res = runGramSchmidt(cols, gsMode);
        const resAlt = runGramSchmidt(cols, gsMode === 'MGS' ? 'CGS' : 'MGS');

        container.querySelector('#m5ALbl').innerHTML = `Editable $A$ (${m}×${n})`;
        container.querySelector('#m5QLbl').innerHTML = `$Q$ (${m}×${n})`;
        container.querySelector('#m5RLbl').innerHTML = `$R$ (${n}×${n}, ${gsMode})`;
        container.querySelector('#m5QOut').innerHTML = renderMatHTML(res.Qmat, 3);
        container.querySelector('#m5ROut').innerHTML = renderMatHTML(res.R, 3);
        renderMath(container.querySelector('#m5ALbl'));
        renderMath(container.querySelector('#m5QLbl'));
        renderMath(container.querySelector('#m5RLbl'));

        const badgesEl = container.querySelector('#m5QrBadges');
        if (res.hasDep) {
          badgesEl.innerHTML = `
            <span class="m5-badge bad">✗ Linearly Dependent Columns ($r_{ii} = 0$)</span>
            <span class="m5-badge warn">Reconstruction $\\|QR - A\\|_{\\max} = ${res.reconErr.toExponential(1)}$</span>
          `;
        } else {
          const orthoClass = res.orthoErr < 1e-6 ? 'ok' : (res.orthoErr < 1e-2 ? 'warn' : 'bad');
          const orthoIcon = orthoClass === 'ok' ? '✓' : '⚠';
          const altLabel = gsMode === 'MGS' ? 'CGS' : 'MGS';
          badgesEl.innerHTML = `
            <span class="m5-badge ${orthoClass}">${orthoIcon} ${gsMode} Orthogonality $\\|Q^TQ - I_${n}\\|_{\\max} = ${res.orthoErr.toExponential(2)}$</span>
            <span class="m5-badge ok">✓ Reconstruction $\\|QR - A\\|_{\\max} = ${res.reconErr.toExponential(2)}$</span>
            <span class="m5-badge ${resAlt.orthoErr < 1e-6 ? 'ok' : 'bad'}">Compare ${altLabel}: $\\|Q^TQ - I_${n}\\|_{\\max} = ${resAlt.orthoErr.toExponential(2)}$</span>
          `;
        }
        renderMath(badgesEl);

        const insEl = container.querySelector('#m5QrInsight');
        insEl.innerHTML = qrPresets[activeQrKey].note;
        renderMath(insEl);
      }

      container.querySelectorAll('[data-qrp]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-qrp]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          activeQrKey = btn.dataset.qrp;
          currentMat = qrPresets[activeQrKey].mat.map(r => r.slice());
          buildAEditor();
          computeQRDisplay();
        });
      });
      container.querySelectorAll('[data-gsmode]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-gsmode]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          gsMode = btn.dataset.gsmode;
          computeQRDisplay();
        });
      });

      // =====================================================================
      // SECTION E: QUICK CONCEPT CHECK
      // =====================================================================
      const explanations = {
        1: 'Just like Module 4, the additive projection formula ∑ (⟨y, uᵢ⟩ / ‖uᵢ‖²) uᵢ equals the true orthogonal projection onto span{u₁, u₂} ONLY when {u₁, u₂} is orthogonal. That is why we must project onto v₁, v₂ (already orthogonalized), NOT x₁, x₂!',
        2: 'v₃ = x₃ - proj_{W₂}(x₃) = 0 means x₃ equals its own projection onto W₂ = span{x₁, x₂}, so x₃ is linearly dependent on x₁ and x₂. Simply discard v₃; {v₁, v₂} is already an orthogonal basis for the entire span.',
        3: 'Left-multiply A = QR by Qᵀ: QᵀA = (QᵀQ)R = IₙR = R. This avoids re-deriving coefficients and guarantees exact consistency with your normalized Q!',
        4: 'In the thin QR factorization of an m×n matrix A (m ≥ n) with independent columns, Q is m×n with QᵀQ = Iₙ (NOTE: QQᵀ is the m×m projection matrix onto Col(A), NOT Iₘ when m > n!), and R is n×n upper triangular with positive diagonal entries rᵢᵢ = ‖vᵢ‖.'
      };

      container.querySelectorAll('.m5-quiz-item').forEach(item => {
        const qNum = item.dataset.q;
        const ans = item.dataset.ans;
        const fb = container.querySelector(`#m5Fb${qNum}`);
        item.querySelectorAll('.m5-opt-btn').forEach(btn => {
          btn.addEventListener('click', () => {
            item.querySelectorAll('.m5-opt-btn').forEach(b => b.classList.remove('correct', 'wrong'));
            if (btn.dataset.opt === ans) {
              btn.classList.add('correct');
              fb.style.color = '#166534';
              fb.textContent = `✓ Correct! ${explanations[qNum]}`;
            } else {
              btn.classList.add('wrong');
              fb.style.color = '#991B1B';
              fb.textContent = `✗ Not quite. Hint: ${explanations[qNum]}`;
            }
          });
        });
      });

      window.addEventListener('resize', drawLab2D);

      // Initial renders
      drawLab2D();
      updateLab3D();
      buildAEditor();
      computeQRDisplay();
    }
  });
})();
