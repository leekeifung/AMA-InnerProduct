/**
 * AMA1751 Linear Algebra — Topic 6: Inner Product
 * Module 6: Least-Squares Problems, Normal Equations & Data Fitting
 */
(function () {
  // Self-contained vector & matrix utilities
  function dot(a, b) {
    let s = 0;
    for (let i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }
  function norm(a) { return Math.sqrt(dot(a, a)); }
  function add(a, b) { return a.map((v, i) => v + b[i]); }
  function sub(a, b) { return a.map((v, i) => v - b[i]); }
  function scale(c, a) { return a.map(v => c * v); }
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
  function matVec(A, x) {
    return A.map(row => dot(row, x));
  }
  function clean(x, d = 3) {
    if (!Number.isFinite(x)) return 'NaN';
    if (Math.abs(x) < 1e-13) return '0';
    if (Math.abs(x) < 1e-3) return x.toExponential(2);
    return String(Number(x.toFixed(d)));
  }
  function fmtVec(v, d = 2) {
    return `(${v.map(x => clean(x, d)).join(', ')})ᵀ`;
  }
  function renderMatHTML(M, d = 3) {
    let html = '<table class="m6-mat-table">';
    for (let i = 0; i < M.length; i++) {
      html += '<tr>';
      for (let j = 0; j < M[0].length; j++) {
        const val = M[i][j];
        const isZero = Math.abs(val) < 1e-13;
        html += `<td class="${isZero ? 'm6-zero-cell' : ''}">${clean(val, d)}</td>`;
      }
      html += '</tr>';
    }
    return html + '</table>';
  }
  function renderColHTML(v, d = 3) {
    return renderMatHTML(v.map(x => [x]), d);
  }

  /**
   * Gaussian elimination with partial pivoting for square system M x = rhs.
   * Uses a tight relative pivot threshold (1e-15) so ill-conditioned systems
   * like Lauchli (eps = 5e-8) demonstrate floating-point inaccuracy while
   * genuinely singular systems (rank-deficient or eps = 1e-8) are flagged.
   */
  function solveLinear(M, rhs) {
    const n = M.length;
    const aug = M.map((row, i) => [...row, rhs[i]]);
    let maxEntry = 0;
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) maxEntry = Math.max(maxEntry, Math.abs(M[i][j]));
    }
    const tol = Math.max(1e-15, maxEntry * 1e-15);

    for (let col = 0; col < n; col++) {
      let piv = col;
      for (let r = col + 1; r < n; r++) {
        if (Math.abs(aug[r][col]) > Math.abs(aug[piv][col])) piv = r;
      }
      if (Math.abs(aug[piv][col]) <= tol) {
        return { ok: false, x: null };
      }
      if (piv !== col) {
        const tmp = aug[col]; aug[col] = aug[piv]; aug[piv] = tmp;
      }
      for (let r = col + 1; r < n; r++) {
        const fac = aug[r][col] / aug[col][col];
        aug[r][col] = 0;
        for (let c = col + 1; c <= n; c++) {
          aug[r][c] -= fac * aug[col][c];
        }
      }
    }

    const x = Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
      let s = aug[i][n];
      for (let j = i + 1; j < n; j++) s -= aug[i][j] * x[j];
      x[i] = s / aug[i][i];
    }
    return { ok: true, x };
  }

  /**
   * Solve least-squares min ||b - A x|| via Normal Equations: (A^T A) x = A^T b.
   */
  function solveNormalEq(A, b) {
    const At = transpose(A);
    const AtA = matMul(At, A);
    const Atb = matVec(At, b);
    const sol = solveLinear(AtA, Atb);
    if (!sol.ok) {
      return { ok: false, AtA, Atb, x: null, resVec: null, resNorm: NaN, sse: NaN, orthoErr: NaN };
    }
    const bHat = matVec(A, sol.x);
    const resVec = sub(b, bHat);
    const resNorm = norm(resVec);
    const Ate = matVec(At, resVec);
    const orthoErr = Math.max(...Ate.map(Math.abs));
    return { ok: true, AtA, Atb, x: sol.x, bHat, resVec, resNorm, sse: resNorm * resNorm, orthoErr };
  }

  /**
   * Solve least-squares min ||b - A x|| via Householder QR: R x = Q^T b.
   * Numerically backward-stable even when A^T A rounds to singular!
   */
  function solveHouseholderQR(A, b) {
    const m = A.length, n = A[0].length;
    const R = A.map(row => row.slice());
    const y = b.slice();

    for (let k = 0; k < Math.min(m, n); k++) {
      let sumSq = 0;
      for (let i = k; i < m; i++) sumSq += R[i][k] * R[i][k];
      const normx = Math.sqrt(sumSq);
      if (normx === 0) continue;

      const sgn = R[k][k] >= 0 ? 1 : -1;
      const alpha = -sgn * normx;
      const u = [];
      u.push(R[k][k] - alpha);
      for (let i = k + 1; i < m; i++) u.push(R[i][k]);

      let uNormSq = 0;
      for (let t = 0; t < u.length; t++) uNormSq += u[t] * u[t];
      if (uNormSq === 0) continue;

      for (let j = k; j < n; j++) {
        let dotUR = 0;
        for (let t = 0; t < u.length; t++) dotUR += u[t] * R[k + t][j];
        const fac = (2 * dotUR) / uNormSq;
        for (let t = 0; t < u.length; t++) R[k + t][j] -= fac * u[t];
      }
      for (let i = k + 1; i < m; i++) R[i][k] = 0;

      let dotUy = 0;
      for (let t = 0; t < u.length; t++) dotUy += u[t] * y[k + t];
      const facy = (2 * dotUy) / uNormSq;
      for (let t = 0; t < u.length; t++) y[k + t] -= facy * u[t];
    }

    // Enforce non-negative diagonal on R (consistent with Module 5 convention)
    const Rtop = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i < m ? R[i][j] : 0)));
    const Qtb = Array.from({ length: n }, (_, i) => (i < m ? y[i] : 0));
    let maxDiag = 0;
    for (let i = 0; i < n; i++) {
      if (Rtop[i][i] < 0) {
        for (let j = i; j < n; j++) Rtop[i][j] = -Rtop[i][j];
        Qtb[i] = -Qtb[i];
      }
      maxDiag = Math.max(maxDiag, Math.abs(Rtop[i][i]));
    }

    const tol = Math.max(1e-14, maxDiag * 1e-14);
    for (let i = 0; i < n; i++) {
      if (Math.abs(Rtop[i][i]) <= tol) {
        return { ok: false, R: Rtop, Qtb, x: null, resVec: null, resNorm: NaN, sse: NaN };
      }
    }

    // Back-substitution on Rtop x = Qtb
    const x = Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
      let s = Qtb[i];
      for (let j = i + 1; j < n; j++) s -= Rtop[i][j] * x[j];
      x[i] = s / Rtop[i][i];
    }

    const bHat = matVec(A, x);
    const resVec = sub(b, bHat);
    const resNorm = norm(resVec);
    return { ok: true, R: Rtop, Qtb, x, bHat, resVec, resNorm, sse: resNorm * resNorm };
  }

  AMA1751.registerModule(6, {
    render(meta) {
      return `
        <style>
          .m6-canvas-wrap {
            position: relative;
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 10px;
            overflow: hidden;
          }
          .m6-canvas {
            display: block;
            width: 100%;
            height: 390px;
            cursor: crosshair;
            touch-action: none;
          }
          .m6-controls {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            margin-bottom: 10px;
            align-items: center;
          }
          .m6-chip {
            padding: 6px 11px;
            border-radius: 999px;
            border: 1px solid var(--border);
            background: #F8FAFC;
            font-size: 0.82rem;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.15s;
          }
          .m6-chip:hover, .m6-chip.active {
            background: var(--polyu-red-soft);
            border-color: #E8B4C4;
            color: var(--polyu-red-dark);
          }
          .m6-toggle-row {
            display: flex;
            flex-wrap: wrap;
            gap: 14px;
            align-items: center;
            margin-bottom: 10px;
            font-size: 0.86rem;
            font-weight: 600;
            color: #334155;
          }
          .m6-toggle-row label {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            cursor: pointer;
          }
          .m6-manual-box {
            background: #FAF5FF;
            border: 1px solid #E9D5FF;
            border-radius: 9px;
            padding: 10px 12px;
            margin-bottom: 10px;
            font-size: 0.85rem;
          }
          .m6-slider-row {
            display: flex;
            align-items: center;
            gap: 10px;
            margin: 5px 0;
          }
          .m6-slider-row input[type="range"] { flex: 1; }
          .m6-slider-lbl {
            min-width: 118px;
            font-weight: 700;
            font-size: 0.84rem;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          }
          .m6-mat-table {
            border-collapse: separate;
            border-spacing: 3px;
            display: inline-table;
            padding: 5px 8px;
            border-left: 2px solid #475569;
            border-right: 2px solid #475569;
            border-radius: 6px;
            background: #FFFFFF;
          }
          .m6-mat-table td {
            min-width: 48px;
            text-align: center;
            padding: 4px 5px;
            border-radius: 4px;
            background: #F8FAFC;
            font-weight: 700;
            font-size: 0.82rem;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          }
          .m6-zero-cell { color: #94A3B8; font-weight: 500 !important; }
          .m6-panel {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 12px 14px;
            font-size: 0.88rem;
            line-height: 1.6;
          }
          .m6-stat-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(115px, 1fr));
            gap: 8px;
            margin: 10px 0;
          }
          .m6-stat-card {
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 8px 10px;
            text-align: center;
          }
          .m6-stat-lbl {
            font-size: 0.72rem;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.04em;
            color: var(--text-muted);
          }
          .m6-stat-val {
            font-size: 1.02rem;
            font-weight: 800;
            color: #0F172A;
            margin-top: 2px;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          }
          .m6-plotly {
            width: 100%;
            height: 420px;
            border: 1px solid var(--border);
            border-radius: 10px;
            background: #FFFFFF;
            overflow: hidden;
          }
          .m6-eq-row {
            display: flex;
            flex-wrap: wrap;
            align-items: center;
            gap: 10px;
            overflow-x: auto;
            padding: 8px 0;
          }
          .m6-eq-block {
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 4px;
          }
          .m6-eq-lbl {
            font-size: 0.74rem;
            font-weight: 700;
            text-transform: uppercase;
            color: var(--text-muted);
          }
          .m6-eq-op {
            font-size: 1.15rem;
            font-weight: 800;
            color: #475569;
            padding-top: 14px;
          }
          .m6-edit-table {
            border-collapse: separate;
            border-spacing: 4px;
            padding: 6px 8px;
            border-left: 2px solid #1D4ED8;
            border-right: 2px solid #1D4ED8;
            border-radius: 6px;
            background: #FFFFFF;
          }
          .m6-edit-input {
            width: 68px;
            padding: 4px 5px;
            border: 1px solid #CBD5E1;
            border-radius: 5px;
            text-align: center;
            font-weight: 700;
            font-size: 0.84rem;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          }
          .m6-compare-card {
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 12px 14px;
          }
          .m6-badge {
            display: inline-block;
            padding: 4px 9px;
            border-radius: 7px;
            font-size: 0.8rem;
            font-weight: 700;
            margin-top: 6px;
          }
          .m6-badge.ok { background: #DCFCE7; color: #166534; border: 1px solid #86EFAC; }
          .m6-badge.warn { background: #FEF3C7; color: #92400E; border: 1px solid #FCD34D; }
          .m6-badge.bad { background: #FEE2E2; color: #991B1B; border: 1px solid #FCA5A5; }
          .m6-quiz-item {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 14px 16px;
            margin-bottom: 12px;
          }
          .m6-opt-group { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 10px; }
          .m6-opt-btn {
            padding: 7px 13px;
            border: 1px solid #CBD5E1;
            border-radius: 8px;
            background: #FFFFFF;
            cursor: pointer;
            font-weight: 600;
            font-size: 0.88rem;
          }
          .m6-opt-btn:hover { background: #F1F5F9; }
          .m6-opt-btn.correct { background: #DCFCE7; border-color: #16A34A; color: #166534; }
          .m6-opt-btn.wrong { background: #FEE2E2; border-color: #DC2626; color: #991B1B; }
          .m6-fb { margin-top: 8px; font-size: 0.88rem; font-weight: 500; }
        </style>

        <!-- HEADER -->
        <div class="module-header">
          <span class="slide-ref">Part 6 · Lecture Notes ${meta.slides}</span>
          <h2>6. Least-Squares Problems, Normal Equations &amp; Data Fitting</h2>
          <p>When an overdetermined system $A\\mathbf{x} = \\mathbf{b}$ has <strong>no exact solution</strong> ($\\mathbf{b} \\notin \\operatorname{Col}(A)$), orthogonal projection finds the closest possible vector $\\hat{\\mathbf{b}} = A\\hat{\\mathbf{x}} \\in \\operatorname{Col}(A)$—turning noisy data into optimal regression curves.</p>
        </div>

        <!-- SECTION A: CORE THEOREMS -->
        <div class="card">
          <h3>The Least-Squares Problem, Normal Equations &amp; $QR$ Shortcut</h3>
          <div class="grid-2">
            <div class="callout definition">
              <strong>Definition &amp; Normal Equations:</strong><br/>
              Given an $m \\times n$ matrix $A$ and $\\mathbf{b} \\in \\mathbb{R}^m$, a <strong>least-squares solution</strong> $\\hat{\\mathbf{x}} \\in \\mathbb{R}^n$ minimizes the error $\\|\\mathbf{b} - A\\mathbf{x}\\|$ (equivalently, the sum of squared errors $\\|\\mathbf{b} - A\\mathbf{x}\\|^2$):
              $$A\\hat{\\mathbf{x}} = \\hat{\\mathbf{b}} = \\operatorname{proj}_{\\operatorname{Col}(A)}(\\mathbf{b})$$
              Since the residual $\\mathbf{e} = \\mathbf{b} - A\\hat{\\mathbf{x}}$ must be orthogonal to every column of $A$ ($A^T\\mathbf{e} = \\mathbf{0}$), $\\hat{\\mathbf{x}}$ satisfies the <strong>Normal Equations</strong>:
              $$A^T A\\,\\hat{\\mathbf{x}} = A^T \\mathbf{b}$$
            </div>

            <div class="callout tip">
              <strong>Uniqueness &amp; Solving via $QR$ Factorization:</strong>
              <ul style="margin:6px 0 0 18px; padding:0; font-size:0.9rem; line-height:1.58;">
                <li><strong>Existence:</strong> A least-squares solution $\\hat{\\mathbf{x}}$ <em>always</em> exists (because $\\hat{\\mathbf{b}} = \\operatorname{proj}_{\\operatorname{Col}(A)}(\\mathbf{b})$ always lies in $\\operatorname{Col}(A)$).</li>
                <li><strong>Uniqueness:</strong> $\\hat{\\mathbf{x}}$ is <em>unique</em> $\\iff$ columns of $A$ are linearly independent $\\iff$ $A^T A$ is invertible:
                  $$\\hat{\\mathbf{x}} = (A^T A)^{-1} A^T \\mathbf{b}$$
                </li>
                <li><strong>$QR$ Shortcut (Numerically Superior!):</strong> If $A = QR$ (thin $QR$, $Q^T Q = I_n$), substitute into $A^T A\\hat{\\mathbf{x}} = A^T\\mathbf{b}$:
                  $$R^T (Q^T Q) R \\hat{\\mathbf{x}} = R^T Q^T \\mathbf{b} \\implies R\\hat{\\mathbf{x}} = Q^T \\mathbf{b}$$
                  Solve instantly by back-substitution without squaring the condition number!</li>
              </ul>
            </div>
          </div>
        </div>

        <!-- SECTION B: INTERACTIVE LAB 1 — 2D REGRESSION & RESIDUAL-SQUARES STUDIO -->
        <div class="card">
          <h3>Interactive Lab 1: 2D Data-Fitting &amp; Visual Residual-Squares Studio</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            <strong>Click</strong> empty space to add a data point, <strong>drag</strong> any point to move it, or <strong>double-click</strong> a point to delete it. Watch how the sum of the amber square areas ($\\text{SSE} = \\sum e_i^2$) responds—and try to beat the optimal line manually!
          </p>

          <div class="m6-controls">
            <strong style="font-size:0.84rem;">Dataset:</strong>
            <button class="m6-chip active" data-ds="linear">Linear Trend</button>
            <button class="m6-chip" data-ds="outlier">High-Leverage Outlier</button>
            <button class="m6-chip" data-ds="quad">Parabolic Curve</button>
            <button class="m6-chip" data-ds="textbook">Textbook $(0,1),(1,1),(2,2),(3,2)$</button>
            <button class="m6-chip" id="m6ClearPt" style="margin-left:auto;">Reset Preset</button>
          </div>

          <div class="m6-controls">
            <strong style="font-size:0.84rem;">Model Basis:</strong>
            <button class="m6-chip active" data-deg="1">Linear ($\\beta_0 + \\beta_1 x$)</button>
            <button class="m6-chip" data-deg="2">Quadratic ($\\beta_0 + \\beta_1 x + \\beta_2 x^2$)</button>
            <button class="m6-chip" data-deg="3">Cubic ($\\beta_0 + \\dots + \\beta_3 x^3$)</button>
          </div>

          <div class="m6-toggle-row">
            <label><input type="checkbox" id="m6ShowSq" checked /> Show Residual Squares (Area $= e_i^2$)</label>
            <label><input type="checkbox" id="m6ManualMode" /> Try to Beat Optimal Line (Manual $\\beta_0, \\beta_1$)</label>
          </div>

          <div id="m6ManualBox" class="m6-manual-box" style="display:none;">
            <div style="font-weight:700; color:#6B21A8; margin-bottom:4px;">Manual Candidate Line (Purple Dashed) vs. Optimal Least-Squares Line (Teal):</div>
            <div class="m6-slider-row">
              <span class="m6-slider-lbl" id="m6ManB0Lbl">β₀ = 0.50</span>
              <input type="range" id="m6ManB0" min="-2" max="8" step="0.05" value="0.5" />
            </div>
            <div class="m6-slider-row">
              <span class="m6-slider-lbl" id="m6ManB1Lbl">β₁ = 0.80</span>
              <input type="range" id="m6ManB1" min="-1.5" max="2.5" step="0.02" value="0.8" />
            </div>
            <div id="m6ManVerdict" style="margin-top:6px; font-weight:700;"></div>
          </div>

          <div class="grid-2" style="align-items:start;">
            <div class="m6-canvas-wrap">
              <canvas id="m6CanvasFit" class="m6-canvas"></canvas>
            </div>

            <div class="m6-panel" id="m6FitStats"></div>
          </div>
        </div>

        <!-- SECTION C: INTERACTIVE LAB 2 — 3D PROJECTION ONTO COL(A) -->
        <div class="card">
          <h3>Interactive Lab 2: Why Normal Equations Work — 3D Geometry of $\\operatorname{Col}(A)$</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            Consider an overdetermined $3\\times 2$ system $A\\mathbf{x} = \\mathbf{b}$. Every possible output $A\\mathbf{x}$ is trapped inside the teal 2D plane $\\operatorname{Col}(A) = \\operatorname{span}\\{\\mathbf{a}_1, \\mathbf{a}_2\\}$. Move $\\mathbf{b}$ or slide a trial vector $\\mathbf{x}_{\\text{trial}}$ to see the <strong>3D right triangle</strong> proving $\\|\\mathbf{b} - A\\mathbf{x}_{\\text{trial}}\\| \\ge \\|\\mathbf{b} - A\\hat{\\mathbf{x}}\\|$!
          </p>

          <div class="m6-controls">
            <strong style="font-size:0.84rem;">$3\\times 2$ System Preset:</strong>
            <button class="m6-chip active" data-p3="clean">Clean Integer Plane ($\\mathbf{b} \\notin \\operatorname{Col}(A)$)</button>
            <button class="m6-chip" data-p3="ortho">Orthogonal Columns $A$</button>
            <button class="m6-chip" data-p3="consistent">Exact Solution ($\\mathbf{b} \\in \\operatorname{Col}(A)$, $\\mathbf{e}=\\mathbf{0}$)</button>
            <button class="m6-chip" id="m6SnapOpt" style="background:#CCFBF1; border-color:#5EEAD4; color:#115E59;">Snap $\\mathbf{x}_{\\text{trial}} \\to \\hat{\\mathbf{x}}$ ✓</button>
          </div>

          <div class="grid-2" style="align-items:start;">
            <div id="m6Plotly3D" class="m6-plotly"></div>

            <div class="m6-panel">
              <div style="font-weight:700; color:#991B1B; margin-bottom:4px;">Target Vector $\\mathbf{b} = (b_1, b_2, b_3)^T$:</div>
              <div class="m6-slider-row">
                <span class="m6-slider-lbl" id="m6B1Lbl">b₁ = 1.0</span>
                <input type="range" id="m6B1" min="-3" max="3" step="0.25" value="1" />
              </div>
              <div class="m6-slider-row">
                <span class="m6-slider-lbl" id="m6B2Lbl">b₂ = 2.0</span>
                <input type="range" id="m6B2" min="-3" max="3" step="0.25" value="2" />
              </div>
              <div class="m6-slider-row">
                <span class="m6-slider-lbl" id="m6B3Lbl">b₃ = 0.0</span>
                <input type="range" id="m6B3" min="-3" max="3" step="0.25" value="0" />
              </div>

              <div style="font-weight:700; color:#7E22CE; margin:10px 0 4px;">Trial Choice $\\mathbf{x}_{\\text{trial}} = (x_1, x_2)^T$ in $\\operatorname{Col}(A)$:</div>
              <div class="m6-slider-row">
                <span class="m6-slider-lbl" id="m6X1Lbl">x₁ = 1.50</span>
                <input type="range" id="m6X1" min="-2.5" max="2.5" step="0.1" value="1.5" />
              </div>
              <div class="m6-slider-row">
                <span class="m6-slider-lbl" id="m6X2Lbl">x₂ = 0.50</span>
                <input type="range" id="m6X2" min="-2.5" max="2.5" step="0.1" value="0.5" />
              </div>

              <div id="m6Readout3D" style="margin-top:10px; border-top:1px solid var(--border); padding-top:10px;"></div>
            </div>
          </div>
        </div>

        <!-- SECTION D: INTERACTIVE LAB 3 — NORMAL EQUATIONS VS QR STABILITY LAB -->
        <div class="card">
          <h3>Interactive Lab 3: Editable Solver — Normal Equations vs. Householder $QR$ Stability</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            Edit any entry of $A$ or $\\mathbf{b}$ below to compare <strong>Normal Equations ($A^T A\\hat{\\mathbf{x}} = A^T\\mathbf{b}$)</strong> against <strong>Householder $QR$ ($R\\hat{\\mathbf{x}} = Q^T\\mathbf{b}$)</strong>. Forming $A^T A$ <em>squares the condition number</em> ($\\kappa(A^T A) = \\kappa(A)^2$): watch what happens on the Lauchli presets when $\\varepsilon = 5\\times 10^{-8}$ and $\\varepsilon = 10^{-8}$!
          </p>

          <div class="m6-controls">
            <strong style="font-size:0.84rem;">Presets:</strong>
            <button class="m6-chip active" data-lsp="exam4x2">$4\\times 2$ Textbook Line Fit</button>
            <button class="m6-chip" data-lsp="dep3x2">$3\\times 2$ Rank-Deficient (Infinitely Many Solutions)</button>
            <button class="m6-chip" data-lsp="lauchliMild">$3\\times 2$ Ill-Conditioned ($\\varepsilon = 5\\times 10^{-8}$)</button>
            <button class="m6-chip" data-lsp="lauchliHard">$3\\times 2$ Lauchli Breakdown ($\\varepsilon = 10^{-8}$)</button>
          </div>

          <div class="m6-eq-row" style="background:#F8FAFC; border:1px solid var(--border); border-radius:10px; padding:12px;">
            <div class="m6-eq-block">
              <div class="m6-eq-lbl" id="m6EditALbl">Editable $A$</div>
              <div id="m6EditA"></div>
            </div>
            <div class="m6-eq-op">$\\mathbf{x} \\approx$</div>
            <div class="m6-eq-block">
              <div class="m6-eq-lbl" id="m6EditBLbl">Editable $\\mathbf{b}$</div>
              <div id="m6EditB"></div>
            </div>
          </div>

          <div class="grid-2" style="margin-top:12px; align-items:start;">
            <div class="m6-compare-card" id="m6CardNE"></div>
            <div class="m6-compare-card" id="m6CardQR"></div>
          </div>

          <div class="callout tip" id="m6LsInsight" style="margin-top:12px;"></div>

          <details class="code-details">
            <summary>Python (SymPy &amp; NumPy) Companion Code — Exact Normal Equations &amp; Stable Least-Squares</summary>
            <pre><code>import sympy as sp
import numpy as np

# 1) Exact symbolic Least-Squares with SymPy (textbook 4x2 fit)
A = sp.Matrix([[1, 0], [1, 1], [1, 2], [1, 3]])
b = sp.Matrix([1, 1, 2, 2])

# Method A: Normal Equations (exact rational arithmetic)
x_hat = (A.T * A).LUsolve(A.T * b)
e = b - A * x_hat
print("Exact x_hat =", x_hat.T)               # Matrix([[9/10, 2/5]])
print("Verify A^T e == 0:", (A.T * e).T)      # Matrix([[0, 0]])
print("Exact SSE =", (e.T * e)[0])            # 1/5

# Method B: QR Least-Squares (exact)
Q, R = A.QRdecomposition()
x_qr = R.LUsolve(Q.T * b)

# 2) Numerical Least-Squares on Lauchli (eps = 1e-8) with NumPy
eps = 1e-8
A_l = np.array([[1.0, 1.0], [eps, 0.0], [0.0, eps]])
b_l = np.array([3.0, eps, 2 * eps])  # Exact solution is [1, 2]

# Normal Equations FAILS because 1 + eps^2 == 1.0 in float64!
try:
    x_ne = np.linalg.solve(A_l.T @ A_l, A_l.T @ b_l)
    print("Normal Eq:", x_ne)
except np.linalg.LinAlgError as err:
    print("Normal Eq FAILED:", err)

# QR / SVD (np.linalg.lstsq) stays exact to 15 digits!
Q_l, R_l = np.linalg.qr(A_l, mode='reduced')
x_qr_np = np.linalg.solve(R_l, Q_l.T @ b_l)
print("Householder QR:", x_qr_np)             # [1. 2.]</code></pre>
          </details>
        </div>

        <!-- SECTION E: QUICK CONCEPT CHECK -->
        <div class="card">
          <h3>Quick Concept Check (Exam Pitfall Radar)</h3>

          <div class="m6-quiz-item" data-q="1" data-ans="C">
            <strong>Q1.</strong> Suppose $A$ is a $6 \\times 3$ matrix whose columns are <em>linearly dependent</em> (so $\\operatorname{rank}(A) = 2$), and $\\mathbf{b} \\in \\mathbb{R}^6$. How many least-squares solutions $\\hat{\\mathbf{x}}$ does $A\\mathbf{x} = \\mathbf{b}$ have?
            <div class="m6-opt-group">
              <button class="m6-opt-btn" data-opt="A">Zero — if $A^T A$ is singular, no least-squares solution exists.</button>
              <button class="m6-opt-btn" data-opt="B">Exactly one, given by $\\hat{\\mathbf{x}} = (A^T A)^{-1} A^T \\mathbf{b}$.</button>
              <button class="m6-opt-btn" data-opt="C">Infinitely many — $\\hat{\\mathbf{b}} = \\operatorname{proj}_{\\operatorname{Col}(A)}(\\mathbf{b})$ is still unique, and any $\\hat{\\mathbf{x}}$ with $A\\hat{\\mathbf{x}} = \\hat{\\mathbf{b}}$ differs by a vector in $\\operatorname{Null}(A)$.</button>
              <button class="m6-opt-btn" data-opt="D">Either zero or one, depending on whether $\\mathbf{b} \\in \\operatorname{Col}(A)$.</button>
            </div>
            <div class="m6-fb" id="m6Fb1"></div>
          </div>

          <div class="m6-quiz-item" data-q="2" data-ans="B">
            <strong>Q2.</strong> If $\\hat{\\mathbf{x}}$ is a least-squares solution of $A\\mathbf{x} = \\mathbf{b}$ with residual $\\mathbf{e} = \\mathbf{b} - A\\hat{\\mathbf{x}}$, which subspace is guaranteed to contain $\\mathbf{e}$?
            <div class="m6-opt-group">
              <button class="m6-opt-btn" data-opt="A">$\\operatorname{Col}(A)$</button>
              <button class="m6-opt-btn" data-opt="B">$\\operatorname{Col}(A)^\\perp = \\operatorname{Null}(A^T)$ (so $A^T \\mathbf{e} = \\mathbf{0}$)</button>
              <button class="m6-opt-btn" data-opt="C">$\\operatorname{Null}(A)$ (so $A\\mathbf{e} = \\mathbf{0}$)</button>
              <button class="m6-opt-btn" data-opt="D">$\\operatorname{Row}(A)$</button>
            </div>
            <div class="m6-fb" id="m6Fb2"></div>
          </div>

          <div class="m6-quiz-item" data-q="3" data-ans="D">
            <strong>Q3.</strong> You want to fit a parabola $y = \\beta_0 + \\beta_1 t + \\beta_2 t^2$ to $5$ data points $(t_i, y_i)$ using least squares. What is the shape of the design matrix $A$, and why is this still a <em>linear</em> least-squares problem?
            <div class="m6-opt-group">
              <button class="m6-opt-btn" data-opt="A">$3 \\times 5$; because $t^2$ is a quadratic function.</button>
              <button class="m6-opt-btn" data-opt="B">$5 \\times 2$; we only need $\\beta_1$ and $\\beta_2$.</button>
              <button class="m6-opt-btn" data-opt="C">It is actually a <em>nonlinear</em> least-squares problem and cannot be solved with Normal Equations.</button>
              <button class="m6-opt-btn" data-opt="D">$5 \\times 3$ (rows $[1, t_i, t_i^2]$); it is linear in the unknown parameter vector $\\boldsymbol{\\beta} = (\\beta_0, \\beta_1, \\beta_2)^T$.</button>
            </div>
            <div class="m6-fb" id="m6Fb3"></div>
          </div>

          <div class="m6-quiz-item" data-q="4" data-ans="A">
            <strong>Q4.</strong> Suppose $A$ is $m \\times n$ with thin $QR$ factorization $A = QR$ ($Q^T Q = I_n$, $R$ invertible upper triangular). Which equation yields the exact least-squares solution $\\hat{\\mathbf{x}}$?
            <div class="m6-opt-group">
              <button class="m6-opt-btn" data-opt="A">$R\\hat{\\mathbf{x}} = Q^T \\mathbf{b}$</button>
              <button class="m6-opt-btn" data-opt="B">$Q\\hat{\\mathbf{x}} = R^T \\mathbf{b}$</button>
              <button class="m6-opt-btn" data-opt="C">$R^T\\hat{\\mathbf{x}} = Q \\mathbf{b}$</button>
              <button class="m6-opt-btn" data-opt="D">$Q^T R\\hat{\\mathbf{x}} = \\mathbf{b}$</button>
            </div>
            <div class="m6-fb" id="m6Fb4"></div>
          </div>
        </div>
      `;
    },

    mount(container, { renderMath }) {
      // =====================================================================
      // LAB 1: 2D DATA-FITTING & RESIDUAL-SQUARES STUDIO
      // =====================================================================
      const canvas = container.querySelector('#m6CanvasFit');
      const ctx = canvas.getContext('2d');
      const datasets = {
        linear:   [[1, 1.6], [2, 2.3], [3, 3.6], [5, 4.7], [6, 6.0], [7, 6.7]],
        outlier:  [[1, 1.5], [2, 2.3], [3, 3.4], [4.5, 4.5], [6, 5.8], [7.4, 1.4]],
        quad:     [[0.6, 6.2], [1.6, 3.5], [2.6, 1.8], [4.0, 1.2], [5.4, 2.0], [6.4, 3.8], [7.4, 6.3]],
        textbook: [[0, 1], [1, 1], [2, 2], [3, 2]]
      };

      let activeDs = 'linear';
      let pts = datasets.linear.map(p => p.slice());
      let degree = 1;
      let showSq = true;
      let manualMode = false;
      let manB0 = 0.5, manB1 = 0.8;
      let dragIdx = -1;

      const xMin = -0.8, xMax = 8.4, yMin = -0.8, yMax = 8.4;

      function fitLayout() {
        const r = canvas.getBoundingClientRect();
        const w = r.width || 460, h = r.height || 390;
        const pad = 34;
        const span = xMax - xMin;
        // Uniform scale so residual squares are true geometric squares!
        const sc = Math.min((w - 2 * pad) / span, (h - 2 * pad) / span);
        const ox = (w - span * sc) / 2 - xMin * sc;
        const oy = h - ((h - span * sc) / 2 - yMin * sc);
        return { w, h, sc, ox, oy };
      }
      function toS(x, y, L) { return [L.ox + x * L.sc, L.oy - y * L.sc]; }
      function toW(sx, sy, L) {
        return [
          Math.max(xMin + 0.2, Math.min(xMax - 0.2, (sx - L.ox) / L.sc)),
          Math.max(yMin + 0.2, Math.min(yMax - 0.2, (L.oy - sy) / L.sc))
        ];
      }

      function buildDesignMatrix(points, deg) {
        const A = points.map(([x]) => {
          const row = [];
          for (let k = 0; k <= deg; k++) row.push(Math.pow(x, k));
          return row;
        });
        const b = points.map(([, y]) => y);
        return { A, b };
      }

      function evalPoly(beta, x) {
        let s = 0;
        for (let k = 0; k < beta.length; k++) s += beta[k] * Math.pow(x, k);
        return s;
      }

      function drawFitLab() {
        if (!canvas.isConnected) return;
        const dpr = window.devicePixelRatio || 1;
        const L = fitLayout();
        canvas.width = Math.round(L.w * dpr);
        canvas.height = Math.round(L.h * dpr);
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.clearRect(0, 0, L.w, L.h);

        // Background grid
        ctx.strokeStyle = '#F1F5F9'; ctx.lineWidth = 1;
        for (let g = 0; g <= 8; g++) {
          const [gx] = toS(g, 0, L), [, gy] = toS(0, g, L);
          ctx.beginPath(); ctx.moveTo(gx, 0); ctx.lineTo(gx, L.h); ctx.stroke();
          ctx.beginPath(); ctx.moveTo(0, gy); ctx.lineTo(L.w, gy); ctx.stroke();
        }
        // Axes
        const [ax0, ay0] = toS(0, 0, L);
        ctx.strokeStyle = '#94A3B8'; ctx.lineWidth = 1.4;
        ctx.beginPath(); ctx.moveTo(0, ay0); ctx.lineTo(L.w, ay0); ctx.stroke();
        ctx.beginPath(); ctx.moveTo(ax0, 0); ctx.lineTo(ax0, L.h); ctx.stroke();

        // Axis tick labels
        ctx.fillStyle = '#64748B'; ctx.font = '11px sans-serif';
        for (let g = 0; g <= 8; g++) {
          const [gx] = toS(g, 0, L), [, gy] = toS(0, g, L);
          ctx.fillText(String(g), gx - 3, ay0 + 14);
          if (g > 0) ctx.fillText(String(g), ax0 - 15, gy + 4);
        }

        const effDeg = Math.min(degree, Math.max(0, pts.length - 1));
        const { A, b } = buildDesignMatrix(pts, effDeg);
        const solQR = solveHouseholderQR(A, b);
        const solNE = solveNormalEq(A, b);
        const beta = solQR.ok ? solQR.x : null;

        // Draw residual squares (either for optimal fit, or for manual line if manualMode is on and deg===1)
        const activeManual = manualMode && degree === 1;
        if (beta && showSq) {
          pts.forEach(([x, y]) => {
            const yhat = activeManual ? (manB0 + manB1 * x) : evalPoly(beta, x);
            const [sx, sy] = toS(x, y, L);
            const [, syh] = toS(x, yhat, L);
            const side = Math.abs(sy - syh);
            const top = Math.min(sy, syh);
            ctx.save();
            ctx.fillStyle = activeManual ? 'rgba(147, 51, 234, 0.14)' : 'rgba(217, 119, 6, 0.16)';
            ctx.strokeStyle = activeManual ? 'rgba(147, 51, 234, 0.55)' : 'rgba(217, 119, 6, 0.65)';
            ctx.lineWidth = 1.2;
            ctx.fillRect(sx, top, side, side);
            ctx.strokeRect(sx, top, side, side);
            ctx.restore();
          });
        }

        // Draw optimal least-squares curve
        if (beta) {
          ctx.save();
          ctx.strokeStyle = '#0F766E'; ctx.lineWidth = 3;
          ctx.beginPath();
          const steps = 120;
          for (let i = 0; i <= steps; i++) {
            const x = xMin + (xMax - xMin) * (i / steps);
            const y = evalPoly(beta, x);
            const [sx, sy] = toS(x, y, L);
            if (i === 0) ctx.moveTo(sx, sy); else ctx.lineTo(sx, sy);
          }
          ctx.stroke();
          ctx.restore();

          // Vertical residual segments for optimal fit
          pts.forEach(([x, y]) => {
            const yhat = evalPoly(beta, x);
            const [sx, sy] = toS(x, y, L);
            const [, syh] = toS(x, yhat, L);
            ctx.save();
            ctx.strokeStyle = activeManual ? 'rgba(15, 118, 110, 0.35)' : '#D97706';
            ctx.lineWidth = 2;
            ctx.setLineDash([4, 3]);
            ctx.beginPath(); ctx.moveTo(sx, sy); ctx.lineTo(sx, syh); ctx.stroke();
            ctx.restore();
          });
        }

        // Draw manual candidate line if enabled
        if (activeManual) {
          const [sx1, sy1] = toS(xMin, manB0 + manB1 * xMin, L);
          const [sx2, sy2] = toS(xMax, manB0 + manB1 * xMax, L);
          ctx.save();
          ctx.strokeStyle = '#7E22CE'; ctx.lineWidth = 2.6; ctx.setLineDash([7, 5]);
          ctx.beginPath(); ctx.moveTo(sx1, sy1); ctx.lineTo(sx2, sy2); ctx.stroke();
          ctx.restore();
        }

        // Draw data points
        pts.forEach(([x, y], idx) => {
          const [sx, sy] = toS(x, y, L);
          ctx.beginPath(); ctx.arc(sx, sy, 6.2, 0, Math.PI * 2);
          ctx.fillStyle = idx === dragIdx ? '#DC2626' : '#1D4ED8';
          ctx.fill();
          ctx.lineWidth = 2; ctx.strokeStyle = '#FFFFFF'; ctx.stroke();
        });

        // Hint tag inside canvas
        ctx.fillStyle = 'rgba(15, 23, 42, 0.58)'; ctx.font = '11px sans-serif';
        ctx.fillText('Click to add · Drag to move · Double-click point to delete', 10, 16);

        // Update manual mode verdict
        if (activeManual && beta) {
          let manSSE = 0;
          pts.forEach(([x, y]) => {
            const e = y - (manB0 + manB1 * x);
            manSSE += e * e;
          });
          const optSSE = solQR.sse;
          const diff = manSSE - optSSE;
          const vEl = container.querySelector('#m6ManVerdict');
          if (diff < 0.015) {
            vEl.style.color = '#15803D';
            vEl.innerHTML = `Your SSE = ${clean(manSSE, 3)} vs. Optimal SSE = ${clean(optSSE, 3)} — 🎯 Spot on the Normal Equations minimum!`;
          } else {
            vEl.style.color = '#6B21A8';
            vEl.innerHTML = `Your SSE = ${clean(manSSE, 3)} vs. Optimal SSE = ${clean(optSSE, 3)} (+${clean(diff, 3)} higher than Normal Equations!)`;
          }
        }

        // Update algebraic panel
        const statsEl = container.querySelector('#m6FitStats');
        if (!beta) {
          statsEl.innerHTML = `<strong style="color:#991B1B;">Rank-Deficient Design Matrix!</strong><br/>Need at least ${degree + 1} points with distinct $x$-coordinates to fit a degree-${degree} polynomial uniquely.`;
          renderMath(statsEl);
          return;
        }

        const yMean = b.reduce((s, v) => s + v, 0) / b.length;
        const sst = b.reduce((s, v) => s + (v - yMean) * (v - yMean), 0);
        const r2 = sst > 1e-12 ? Math.max(0, 1 - solQR.sse / sst) : 1;
        const rmse = Math.sqrt(solQR.sse / pts.length);
        const At = transpose(A);
        const Ate = matVec(At, solQR.resVec);
        const orthoMax = Math.max(...Ate.map(Math.abs));

        const terms = beta.map((c, k) => {
          const val = clean(c, 3);
          if (k === 0) return val;
          const sign = c >= 0 ? ' + ' : ' - ';
          const absVal = clean(Math.abs(c), 3);
          return `${sign}${absVal}x${k > 1 ? `^${k}` : ''}`;
        }).join('');

        statsEl.innerHTML = `
          <div><strong>Fitted Model ($n = ${pts.length}$ points, $A$ is $${pts.length}\\times ${effDeg + 1}$):</strong></div>
          <div style="margin:4px 0 8px; padding:7px 10px; background:#F0FDFA; border-left:3px solid #0F766E; border-radius:6px; font-weight:700; color:#115E59;">
            $\\hat{y} = ${terms}$
          </div>
          <div class="m6-stat-grid">
            <div class="m6-stat-card"><div class="m6-stat-lbl">SSE $\\|\\mathbf{e}\\|^2$</div><div class="m6-stat-val">${clean(solQR.sse, 3)}</div></div>
            <div class="m6-stat-card"><div class="m6-stat-lbl">RMSE</div><div class="m6-stat-val">${clean(rmse, 3)}</div></div>
            <div class="m6-stat-card"><div class="m6-stat-lbl">$R^2$ Fit</div><div class="m6-stat-val">${clean(r2, 3)}</div></div>
            <div class="m6-stat-card"><div class="m6-stat-lbl">$\\|A^T\\mathbf{e}\\|_\\infty$</div><div class="m6-stat-val" style="color:#15803D;">${orthoMax < 1e-10 ? '0 ✓' : orthoMax.toExponential(1)}</div></div>
          </div>
          <div class="m6-eq-row">
            <div class="m6-eq-block"><div class="m6-eq-lbl">$A^T A$</div>${renderMatHTML(solNE.AtA, 2)}</div>
            <div class="m6-eq-block"><div class="m6-eq-lbl">$\\hat{\\boldsymbol{\\beta}}$</div>${renderColHTML(beta, 3)}</div>
            <div class="m6-eq-op">=</div>
            <div class="m6-eq-block"><div class="m6-eq-lbl">$A^T \\mathbf{b}$</div>${renderColHTML(solNE.Atb, 2)}</div>
          </div>
          <div style="font-size:0.8rem; color:var(--text-muted); margin-top:4px;">
            Residual orthogonality verified: $\\sum e_i = ${clean(Ate[0], 3)}$ and $\\sum x_i e_i = ${Ate.length > 1 ? clean(Ate[1], 3) : '0'}$ ✓
          </div>
        `;
        renderMath(statsEl);
      }

      function canvasPtr(e) {
        const r = canvas.getBoundingClientRect();
        const cx = e.touches ? e.touches[0].clientX : e.clientX;
        const cy = e.touches ? e.touches[0].clientY : e.clientY;
        return [cx - r.left, cy - r.top];
      }
      function findHit(px, py, L) {
        for (let i = 0; i < pts.length; i++) {
          const [sx, sy] = toS(pts[i][0], pts[i][1], L);
          if (Math.hypot(px - sx, py - sy) < 14) return i;
        }
        return -1;
      }

      canvas.addEventListener('mousedown', e => {
        const L = fitLayout(), [px, py] = canvasPtr(e);
        const hit = findHit(px, py, L);
        if (hit >= 0) {
          dragIdx = hit;
        } else {
          const [wx, wy] = toW(px, py, L);
          pts.push([Math.round(wx * 10) / 10, Math.round(wy * 10) / 10]);
          dragIdx = pts.length - 1;
          drawFitLab();
        }
      });
      canvas.addEventListener('touchstart', e => {
        const L = fitLayout(), [px, py] = canvasPtr(e);
        const hit = findHit(px, py, L);
        if (hit >= 0) {
          dragIdx = hit; e.preventDefault();
        } else {
          const [wx, wy] = toW(px, py, L);
          pts.push([Math.round(wx * 10) / 10, Math.round(wy * 10) / 10]);
          dragIdx = pts.length - 1;
          drawFitLab();
        }
      }, { passive: false });
      canvas.addEventListener('dblclick', e => {
        if (pts.length <= 2) return;
        const L = fitLayout(), [px, py] = canvasPtr(e);
        const hit = findHit(px, py, L);
        if (hit >= 0) {
          pts.splice(hit, 1);
          dragIdx = -1;
          drawFitLab();
        }
      });
      function onMoveFit(e) {
        if (dragIdx < 0 || !canvas.isConnected) return;
        if (e.touches) e.preventDefault();
        const L = fitLayout(), [px, py] = canvasPtr(e);
        const [wx, wy] = toW(px, py, L);
        pts[dragIdx] = [Math.round(wx * 10) / 10, Math.round(wy * 10) / 10];
        drawFitLab();
      }
      window.addEventListener('mousemove', onMoveFit);
      window.addEventListener('touchmove', onMoveFit, { passive: false });
      window.addEventListener('mouseup', () => { dragIdx = -1; });
      window.addEventListener('touchend', () => { dragIdx = -1; });

      container.querySelectorAll('[data-ds]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-ds]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          activeDs = btn.dataset.ds;
          pts = datasets[activeDs].map(p => p.slice());
          if (activeDs === 'quad') {
            degree = 2;
            container.querySelectorAll('[data-deg]').forEach(b => b.classList.toggle('active', b.dataset.deg === '2'));
            manualMode = false;
            container.querySelector('#m6ManualMode').checked = false;
            container.querySelector('#m6ManualBox').style.display = 'none';
          }
          drawFitLab();
        });
      });
      container.querySelector('#m6ClearPt').addEventListener('click', () => {
        pts = datasets[activeDs].map(p => p.slice());
        drawFitLab();
      });
      container.querySelectorAll('[data-deg]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-deg]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          degree = parseInt(btn.dataset.deg, 10);
          if (degree !== 1 && manualMode) {
            manualMode = false;
            container.querySelector('#m6ManualMode').checked = false;
            container.querySelector('#m6ManualBox').style.display = 'none';
          }
          drawFitLab();
        });
      });
      container.querySelector('#m6ShowSq').addEventListener('change', e => {
        showSq = e.target.checked;
        drawFitLab();
      });
      container.querySelector('#m6ManualMode').addEventListener('change', e => {
        manualMode = e.target.checked;
        if (manualMode && degree !== 1) {
          degree = 1;
          container.querySelectorAll('[data-deg]').forEach(b => b.classList.toggle('active', b.dataset.deg === '1'));
        }
        container.querySelector('#m6ManualBox').style.display = manualMode ? 'block' : 'none';
        drawFitLab();
      });
      container.querySelector('#m6ManB0').addEventListener('input', e => {
        manB0 = parseFloat(e.target.value);
        container.querySelector('#m6ManB0Lbl').textContent = `β₀ = ${manB0.toFixed(2)}`;
        drawFitLab();
      });
      container.querySelector('#m6ManB1').addEventListener('input', e => {
        manB1 = parseFloat(e.target.value);
        container.querySelector('#m6ManB1Lbl').textContent = `β₁ = ${manB1.toFixed(2)}`;
        drawFitLab();
      });

      // =====================================================================
      // LAB 2: 3D PROJECTION ONTO COL(A) (PLOTLY 3D)
      // =====================================================================
      const p3Presets = {
        clean:      { A: [[1, 0], [0, 1], [1, 1]], b: [1, 2, 0],   xt: [1.5, 0.5] },
        ortho:      { A: [[1, 1], [1, -1], [0, 1]], b: [2, 0, 2.5], xt: [2.0, -0.5] },
        consistent: { A: [[1, 0], [0, 1], [1, 1]], b: [1, 1.5, 2.5], xt: [-0.5, 2.0] }
      };
      let activeP3 = 'clean';
      let A3 = p3Presets.clean.A.map(r => r.slice());
      let b3 = p3Presets.clean.b.slice();
      let xt3 = p3Presets.clean.xt.slice();
      const p3El = container.querySelector('#m6Plotly3D');

      function sync3DSliders() {
        ['B1', 'B2', 'B3'].forEach((id, idx) => {
          container.querySelector(`#m6${id}`).value = b3[idx];
          container.querySelector(`#m6${id}Lbl`).textContent = `b${'₁₂₃'[idx]} = ${b3[idx].toFixed(2)}`;
        });
        ['X1', 'X2'].forEach((id, idx) => {
          container.querySelector(`#m6${id}`).value = xt3[idx];
          container.querySelector(`#m6${id}Lbl`).textContent = `x${'₁₂'[idx]} = ${xt3[idx].toFixed(2)}`;
        });
      }

      function update3DLab() {
        const sol = solveNormalEq(A3, b3);
        const xHat = sol.x;
        const bHat = sol.bHat;
        const eVec = sol.resVec;
        const axTrial = matVec(A3, xt3);
        const trialErrVec = sub(b3, axTrial);
        const planeDiff = sub(axTrial, bHat);

        const errOptSq = dot(eVec, eVec);
        const errTrialSq = dot(trialErrVec, trialErrVec);
        const diffSq = dot(planeDiff, planeDiff);

        const ro = container.querySelector('#m6Readout3D');
        ro.innerHTML = `
          <div><strong>Normal Equations Solution:</strong> $\\hat{\\mathbf{x}} = ${fmtVec(xHat, 2)}$, &nbsp;$\\hat{\\mathbf{b}} = A\\hat{\\mathbf{x}} = ${fmtVec(bHat, 2)}$</div>
          <div><strong>Optimal Residual:</strong> $\\mathbf{e} = \\mathbf{b} - \\hat{\\mathbf{b}} = ${fmtVec(eVec, 2)}$, &nbsp;$A^T\\mathbf{e} = ${fmtVec(matVec(transpose(A3), eVec), 3)}$ ✓</div>
          <div style="margin-top:7px; padding:8px 10px; background:#FFFFFF; border:1px solid var(--border); border-radius:8px; font-size:0.84rem;">
            <strong>Pythagorean Optimality Proof:</strong><br/>
            $\\|\\mathbf{b} - A\\mathbf{x}_{\\text{trial}}\\|^2 = \\|\\mathbf{e}\\|^2 + \\|A\\hat{\\mathbf{x}} - A\\mathbf{x}_{\\text{trial}}\\|^2$<br/>
            <span style="color:#7E22CE; font-weight:800;">${clean(errTrialSq, 3)}</span>
            = <span style="color:#D97706; font-weight:800;">${clean(errOptSq, 3)}</span>
            + <span style="color:#0F766E; font-weight:800;">${clean(diffSq, 3)}</span>
            ${diffSq < 1e-4 ? ' &nbsp;<strong style="color:#15803D;">(Minimum Reached! ✓)</strong>' : ''}
          </div>
        `;
        renderMath(ro);

        if (typeof Plotly === 'undefined') {
          p3El.innerHTML = '<div style="padding:20px; color:var(--text-muted);">3D plot requires Plotly.js.</div>';
          return;
        }

        const a1 = [A3[0][0], A3[1][0], A3[2][0]];
        const a2 = [A3[0][1], A3[1][1], A3[2][1]];

        // Build plane surface for Col(A) = span{a1, a2}
        const gv = [-2.6, 0, 2.6];
        const sx = [], sy = [], sz = [];
        for (const u of gv) {
          const rx = [], ry = [], rz = [];
          for (const v of gv) {
            const pt = add(scale(u, a1), scale(v, a2));
            rx.push(pt[0]); ry.push(pt[1]); rz.push(pt[2]);
          }
          sx.push(rx); sy.push(ry); sz.push(rz);
        }

        const traces = [
          {
            type: 'surface', name: 'Col(A) plane', showscale: false, opacity: 0.28,
            x: sx, y: sy, z: sz,
            colorscale: [[0, '#99F6E4'], [1, '#0F766E']], hoverinfo: 'skip'
          },
          {
            type: 'scatter3d', mode: 'lines+markers+text', name: 'a₁, a₂ (cols of A)',
            x: [0, a1[0], null, 0, a2[0]], y: [0, a1[1], null, 0, a2[1]], z: [0, a1[2], null, 0, a2[2]],
            text: ['', 'a₁', '', '', 'a₂'], textposition: 'top center',
            textfont: { color: '#1D4ED8', size: 12 },
            line: { color: '#1D4ED8', width: 4 }, marker: { size: 3.5, color: '#1D4ED8' }
          },
          {
            type: 'scatter3d', mode: 'lines+markers+text', name: 'Target b',
            x: [0, b3[0]], y: [0, b3[1]], z: [0, b3[2]],
            text: ['', 'b'], textposition: 'top center',
            textfont: { color: '#DC2626', size: 13 },
            line: { color: '#DC2626', width: 6 }, marker: { size: [2, 6], color: '#DC2626' }
          },
          {
            type: 'scatter3d', mode: 'lines+markers+text', name: 'Optimal b̂ = Ax̂',
            x: [0, bHat[0]], y: [0, bHat[1]], z: [0, bHat[2]],
            text: ['', 'b̂ = Ax̂'], textposition: 'bottom center',
            textfont: { color: '#0F766E', size: 13 },
            line: { color: '#0F766E', width: 6 }, marker: { size: [2, 6], color: '#0F766E' }
          },
          {
            type: 'scatter3d', mode: 'lines+text', name: 'Residual e ⊥ Col(A)',
            x: [bHat[0], b3[0]], y: [bHat[1], b3[1]], z: [bHat[2], b3[2]],
            text: ['', ''],
            line: { color: '#D97706', width: 5, dash: 'dash' }
          },
          {
            type: 'scatter3d', mode: 'markers+text', name: 'Trial Ax',
            x: [axTrial[0]], y: [axTrial[1]], z: [axTrial[2]],
            text: ['Ax_trial'], textposition: 'bottom center',
            textfont: { color: '#7E22CE', size: 12 },
            marker: { size: 6.5, color: '#7E22CE' }
          },
          {
            type: 'scatter3d', mode: 'lines', name: 'Trial error ||b - Ax||',
            x: [axTrial[0], b3[0]], y: [axTrial[1], b3[1]], z: [axTrial[2], b3[2]],
            line: { color: '#7E22CE', width: 4, dash: 'dot' }
          },
          {
            type: 'scatter3d', mode: 'lines', showlegend: false,
            x: [bHat[0], axTrial[0]], y: [bHat[1], axTrial[1]], z: [bHat[2], axTrial[2]],
            line: { color: '#0F766E', width: 3.5 }
          }
        ];

        Plotly.react(p3El, traces, {
          margin: { l: 0, r: 0, t: 6, b: 0 },
          legend: { orientation: 'h', y: 0.99, x: 0.01, font: { size: 10 } },
          scene: {
            aspectmode: 'cube',
            xaxis: { range: [-4, 4], title: 'x₁', zerolinecolor: '#94A3B8' },
            yaxis: { range: [-4, 4], title: 'x₂', zerolinecolor: '#94A3B8' },
            zaxis: { range: [-4, 4], title: 'x₃', zerolinecolor: '#94A3B8' },
            camera: { eye: { x: 1.5, y: 1.45, z: 1.15 } }
          }
        }, { responsive: true, displayModeBar: false });
      }

      ['B1', 'B2', 'B3'].forEach((id, idx) => {
        container.querySelector(`#m6${id}`).addEventListener('input', e => {
          b3[idx] = parseFloat(e.target.value);
          container.querySelector(`#m6${id}Lbl`).textContent = `b${'₁₂₃'[idx]} = ${b3[idx].toFixed(2)}`;
          update3DLab();
        });
      });
      ['X1', 'X2'].forEach((id, idx) => {
        container.querySelector(`#m6${id}`).addEventListener('input', e => {
          xt3[idx] = parseFloat(e.target.value);
          container.querySelector(`#m6${id}Lbl`).textContent = `x${'₁₂'[idx]} = ${xt3[idx].toFixed(2)}`;
          update3DLab();
        });
      });
      container.querySelector('#m6SnapOpt').addEventListener('click', () => {
        const sol = solveNormalEq(A3, b3);
        if (sol.ok) {
          xt3 = sol.x.map(v => Math.max(-2.5, Math.min(2.5, v)));
          sync3DSliders();
          update3DLab();
        }
      });
      container.querySelectorAll('[data-p3]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-p3]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          activeP3 = btn.dataset.p3;
          A3 = p3Presets[activeP3].A.map(r => r.slice());
          b3 = p3Presets[activeP3].b.slice();
          xt3 = p3Presets[activeP3].xt.slice();
          sync3DSliders();
          update3DLab();
        });
      });

      // =====================================================================
      // LAB 3: EDITABLE SOLVER — NORMAL EQUATIONS VS HOUSEHOLDER QR
      // =====================================================================
      const lsPresets = {
        exam4x2: {
          A: [[1, 0], [1, 1], [1, 2], [1, 3]],
          b: [1, 1, 2, 2],
          note: '<strong>Well-Conditioned Textbook System:</strong> Both Normal Equations and Householder $QR$ produce the exact same solution $\\hat{\\mathbf{x}} = (0.9, 0.4)^T$ with minimum error $\\|\\mathbf{b} - A\\hat{\\mathbf{x}}\\| = \\sqrt{0.2} \\approx 0.4472$.'
        },
        dep3x2: {
          A: [[1, 2], [2, 4], [-1, -2]],
          b: [3, 6, 0],
          note: '<strong>Rank-Deficient System ($\\mathbf{a}_2 = 2\\mathbf{a}_1$):</strong> Both $A^T A$ and $R$ are singular! That does <em>not</em> mean no least-squares solution exists—instead, $\\hat{\\mathbf{b}} = (2.5, 5, -2.5)^T$ is still unique, and <strong>every</strong> $(x_1, x_2)^T$ satisfying $x_1 + 2x_2 = 2.5$ is a valid least-squares solution.'
        },
        lauchliMild: {
          A: [[1, 1], [5e-8, 0], [0, 5e-8]],
          b: [3, 5e-8, 1e-7],
          note: '<strong>Squaring the Condition Number ($\\varepsilon = 5\\times 10^{-8}$):</strong> The exact solution is $\\mathbf{x}^* = (1, 2)^T$ with zero residual! Because $1 + \\varepsilon^2 = 1 + 2.5\\times 10^{-15}$ barely fits inside 64-bit IEEE floating point, Normal Equations loses almost all precision ($\\hat{\\mathbf{x}}_{\\text{NE}} \\approx (1.0455, 1.9545)^T$, a ~4.5% error!), whereas Householder $QR$ recovers $(1, 2)^T$ cleanly.'
        },
        lauchliHard: {
          A: [[1, 1], [1e-8, 0], [0, 1e-8]],
          b: [3, 1e-8, 2e-8],
          note: '<strong>Catastrophic Breakdown of Normal Equations ($\\varepsilon = 10^{-8}$):</strong> In 64-bit floating-point arithmetic, $1 + (10^{-8})^2 = 1 + 10^{-16} \\equiv 1.0$, so $A^T A$ rounds to $\\begin{bmatrix}1 & 1\\\\ 1 & 1\\end{bmatrix}$ and becomes <strong>numerically singular</strong>! Meanwhile, Householder $QR$ avoids forming $A^T A$ and solves $R\\hat{\\mathbf{x}} = Q^T\\mathbf{b}$ to exact $(1, 2)^T$!'
        }
      };

      let activeLsp = 'exam4x2';
      let curA = lsPresets.exam4x2.A.map(r => r.slice());
      let curB = lsPresets.exam4x2.b.slice();

      function fmtInputVal(v) {
        if (v === 0) return '0';
        if (Math.abs(v) < 1e-4) return v.toExponential(0);
        return String(v);
      }

      function buildLsEditors() {
        const m = curA.length, n = curA[0].length;
        container.querySelector('#m6EditALbl').innerHTML = `Editable $A$ (${m}×${n})`;
        container.querySelector('#m6EditBLbl').innerHTML = `Editable $\\mathbf{b}$ (${m}×1)`;
        renderMath(container.querySelector('#m6EditALbl'));
        renderMath(container.querySelector('#m6EditBLbl'));

        let htmlA = '<table class="m6-edit-table">';
        for (let i = 0; i < m; i++) {
          htmlA += '<tr>';
          for (let j = 0; j < n; j++) {
            htmlA += `<td><input type="number" step="any" class="m6-edit-input m6-in-a" data-r="${i}" data-c="${j}" value="${fmtInputVal(curA[i][j])}" /></td>`;
          }
          htmlA += '</tr>';
        }
        htmlA += '</table>';
        const edA = container.querySelector('#m6EditA');
        edA.innerHTML = htmlA;

        let htmlB = '<table class="m6-edit-table" style="border-color:#DC2626;">';
        for (let i = 0; i < m; i++) {
          htmlB += `<tr><td><input type="number" step="any" class="m6-edit-input m6-in-b" data-r="${i}" value="${fmtInputVal(curB[i])}" /></td></tr>`;
        }
        htmlB += '</table>';
        const edB = container.querySelector('#m6EditB');
        edB.innerHTML = htmlB;

        edA.querySelectorAll('.m6-in-a').forEach(inp => {
          inp.addEventListener('input', () => {
            const v = parseFloat(inp.value);
            curA[parseInt(inp.dataset.r, 10)][parseInt(inp.dataset.c, 10)] = Number.isFinite(v) ? v : 0;
            updateLsCompare();
          });
        });
        edB.querySelectorAll('.m6-in-b').forEach(inp => {
          inp.addEventListener('input', () => {
            const v = parseFloat(inp.value);
            curB[parseInt(inp.dataset.r, 10)] = Number.isFinite(v) ? v : 0;
            updateLsCompare();
          });
        });
      }

      function updateLsCompare() {
        const resNE = solveNormalEq(curA, curB);
        const resQR = solveHouseholderQR(curA, curB);

        // Card 1: Normal Equations
        const cNE = container.querySelector('#m6CardNE');
        if (!resNE.ok) {
          cNE.innerHTML = `
            <h4 style="margin:0 0 8px; color:#991B1B;">Method 1: Normal Equations ($A^T A\\hat{\\mathbf{x}} = A^T\\mathbf{b}$)</h4>
            <div class="m6-eq-row">
              <div class="m6-eq-block"><div class="m6-eq-lbl">$A^T A$</div>${renderMatHTML(resNE.AtA, 4)}</div>
              <div class="m6-eq-block"><div class="m6-eq-lbl">$A^T \\mathbf{b}$</div>${renderColHTML(resNE.Atb, 4)}</div>
            </div>
            <div class="m6-badge bad">✗ FAILED: $A^T A$ is singular (numerically or structurally)!</div>
          `;
        } else {
          const diffFromQR = resQR.ok ? norm(sub(resNE.x, resQR.x)) : 0;
          const badgeClass = diffFromQR > 1e-3 ? 'warn' : 'ok';
          const badgeTxt = diffFromQR > 1e-3
            ? `⚠ Lost Precision! $\\|\\hat{\\mathbf{x}}_{\\text{NE}} - \\hat{\\mathbf{x}}_{QR}\\| = ${clean(diffFromQR, 4)}$`
            : `✓ Residual $\\|\\mathbf{b} - A\\hat{\\mathbf{x}}\\| = ${clean(resNE.resNorm, 4)}$`;
          cNE.innerHTML = `
            <h4 style="margin:0 0 8px; color:#1E3A8A;">Method 1: Normal Equations ($A^T A\\hat{\\mathbf{x}} = A^T\\mathbf{b}$)</h4>
            <div class="m6-eq-row">
              <div class="m6-eq-block"><div class="m6-eq-lbl">$A^T A$</div>${renderMatHTML(resNE.AtA, 3)}</div>
              <div class="m6-eq-block"><div class="m6-eq-lbl">$\\hat{\\mathbf{x}}_{\\text{NE}}$</div>${renderColHTML(resNE.x, 4)}</div>
              <div class="m6-eq-op">=</div>
              <div class="m6-eq-block"><div class="m6-eq-lbl">$A^T \\mathbf{b}$</div>${renderColHTML(resNE.Atb, 3)}</div>
            </div>
            <div class="m6-badge ${badgeClass}">${badgeTxt}</div>
          `;
        }

        // Card 2: Householder QR
        const cQR = container.querySelector('#m6CardQR');
        if (!resQR.ok) {
          cQR.innerHTML = `
            <h4 style="margin:0 0 8px; color:#0F766E;">Method 2: Householder $QR$ ($R\\hat{\\mathbf{x}} = Q^T\\mathbf{b}$)</h4>
            <div class="m6-eq-row">
              <div class="m6-eq-block"><div class="m6-eq-lbl">$R$</div>${renderMatHTML(resQR.R, 3)}</div>
              <div class="m6-eq-block"><div class="m6-eq-lbl">$Q^T \\mathbf{b}$</div>${renderColHTML(resQR.Qtb, 3)}</div>
            </div>
            <div class="m6-badge warn">⚠ Rank-Deficient $A$ ($r_{ii} = 0$): Infinitely many LS solutions!</div>
          `;
        } else {
          cQR.innerHTML = `
            <h4 style="margin:0 0 8px; color:#0F766E;">Method 2: Householder $QR$ ($R\\hat{\\mathbf{x}} = Q^T\\mathbf{b}$)</h4>
            <div class="m6-eq-row">
              <div class="m6-eq-block"><div class="m6-eq-lbl">$R$</div>${renderMatHTML(resQR.R, 3)}</div>
              <div class="m6-eq-block"><div class="m6-eq-lbl">$\\hat{\\mathbf{x}}_{QR}$</div>${renderColHTML(resQR.x, 4)}</div>
              <div class="m6-eq-op">=</div>
              <div class="m6-eq-block"><div class="m6-eq-lbl">$Q^T \\mathbf{b}$</div>${renderColHTML(resQR.Qtb, 3)}</div>
            </div>
            <div class="m6-badge ok">✓ Stable! Residual $\\|\\mathbf{b} - A\\hat{\\mathbf{x}}_{QR}\\| = ${clean(resQR.resNorm, 4)}$</div>
          `;
        }

        const ins = container.querySelector('#m6LsInsight');
        ins.innerHTML = lsPresets[activeLsp].note;
        renderMath(cNE); renderMath(cQR); renderMath(ins);
      }

      container.querySelectorAll('[data-lsp]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-lsp]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          activeLsp = btn.dataset.lsp;
          curA = lsPresets[activeLsp].A.map(r => r.slice());
          curB = lsPresets[activeLsp].b.slice();
          buildLsEditors();
          updateLsCompare();
        });
      });

      // =====================================================================
      // SECTION E: QUICK CONCEPT CHECK
      // =====================================================================
      const explanations = {
        1: 'A least-squares solution ALWAYS exists because b̂ = proj_{Col(A)}(b) always lies in Col(A). When the columns of A are linearly dependent, Null(A) is nontrivial, so adding any vector z ∈ Null(A) to x̂ gives A(x̂ + z) = Ax̂ = b̂—producing infinitely many least-squares solutions!',
        2: 'Because b̂ = Ax̂ is the orthogonal projection of b onto Col(A), the residual e = b - Ax̂ is perpendicular to every column of A. That is literally the matrix equation Aᵀe = 0, i.e. e ∈ Null(Aᵀ) = Col(A)⊥.',
        3: 'The word "linear" in linear least squares refers to linearity in the UNKNOWN PARAMETERS β₀, β₁, β₂—NOT in the input variable t! Each row of the 5×3 design matrix A is [1, tᵢ, tᵢ²], so Aβ = y is a standard linear system in β.',
        4: 'Substitute A = QR into the Normal Equations AᵀAx̂ = Aᵀb: Rᵀ(QᵀQ)Rx̂ = RᵀQᵀb. Since QᵀQ = Iₙ and Rᵀ is invertible, cancel Rᵀ from both sides to obtain Rx̂ = Qᵀb!'
      };

      container.querySelectorAll('.m6-quiz-item').forEach(item => {
        const qNum = item.dataset.q;
        const ans = item.dataset.ans;
        const fb = container.querySelector(`#m6Fb${qNum}`);
        item.querySelectorAll('.m6-opt-btn').forEach(btn => {
          btn.addEventListener('click', () => {
            item.querySelectorAll('.m6-opt-btn').forEach(b => b.classList.remove('correct', 'wrong'));
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

      window.addEventListener('resize', drawFitLab);

      // Initial renders
      drawFitLab();
      sync3DSliders();
      update3DLab();
      buildLsEditors();
      updateLsCompare();
    }
  });
})();
