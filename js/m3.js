/**
 * AMA1751 Linear Algebra — Topic 6: Inner Product
 * Module 3: Orthogonal/Orthonormal Sets & Orthogonal Matrices
 */
(function () {
  // Self-contained matrix & vector utilities
  function dot(a, b) {
    let s = 0;
    for (let i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }
  function norm(a) { return Math.sqrt(dot(a, a)); }
  function scale(c, a) { return a.map(v => c * v); }
  function add(a, b) { return a.map((v, i) => v + b[i]); }
  function transpose(M) {
    const r = M.length, c = M[0].length;
    const T = Array.from({ length: c }, () => Array(r).fill(0));
    for (let i = 0; i < r; i++)
      for (let j = 0; j < c; j++) T[j][i] = M[i][j];
    return T;
  }
  function matMul(A, B) {
    const r = A.length, k = A[0].length, c = B[0].length;
    const C = Array.from({ length: r }, () => Array(c).fill(0));
    for (let i = 0; i < r; i++)
      for (let j = 0; j < c; j++) {
        let s = 0;
        for (let t = 0; t < k; t++) s += A[i][t] * B[t][j];
        C[i][j] = s;
      }
    return C;
  }
  function matVec2(M, v) {
    return [M[0][0] * v[0] + M[0][1] * v[1], M[1][0] * v[0] + M[1][1] * v[1]];
  }
  function det3(M) {
    return (
      M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) -
      M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
      M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])
    );
  }
  function cleanNum(x, digits = 3) {
    if (Math.abs(x) < 1e-6) return '0';
    const r = Number(x.toFixed(digits));
    return String(r);
  }

  AMA1751.registerModule(3, {
    render(meta) {
      return `
        <style>
          .m3-canvas-wrap {
            position: relative;
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 10px;
            overflow: hidden;
          }
          .m3-canvas {
            display: block;
            width: 100%;
            height: 360px;
            cursor: grab;
            touch-action: none;
          }
          .m3-canvas:active { cursor: grabbing; }
          .m3-canvas-static {
            display: block;
            width: 100%;
            height: 320px;
          }
          .m3-controls {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            margin-bottom: 12px;
            align-items: center;
          }
          .m3-chip {
            padding: 6px 11px;
            border-radius: 999px;
            border: 1px solid var(--border);
            background: #F8FAFC;
            font-size: 0.82rem;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.15s;
          }
          .m3-chip:hover, .m3-chip.active {
            background: var(--polyu-red-soft);
            border-color: #E8B4C4;
            color: var(--polyu-red-dark);
          }
          .m3-stat-grid {
            display: grid;
            grid-template-columns: repeat(2, minmax(0, 1fr));
            gap: 10px;
          }
          .m3-stat {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 10px 12px;
          }
          .m3-stat-label {
            font-size: 0.74rem;
            color: var(--text-muted);
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.04em;
          }
          .m3-stat-val {
            font-size: 0.96rem;
            font-weight: 700;
            margin-top: 3px;
            color: #0F172A;
          }
          .m3-mat-grid {
            display: inline-grid;
            gap: 6px;
            padding: 10px 14px;
            border-left: 2.5px solid #334155;
            border-right: 2.5px solid #334155;
            border-radius: 8px;
            background: #FFFFFF;
          }
          .m3-mat-input {
            width: 68px;
            padding: 6px 6px;
            border: 1px solid #CBD5E1;
            border-radius: 6px;
            text-align: center;
            font-weight: 600;
            font-size: 0.9rem;
          }
          .m3-mat-input:focus {
            outline: 2px solid var(--polyu-red);
            border-color: var(--polyu-red);
          }
          .m3-gram-table {
            border-collapse: separate;
            border-spacing: 4px;
            display: inline-table;
            padding: 6px 10px;
            border-left: 2px solid #475569;
            border-right: 2px solid #475569;
            border-radius: 6px;
            background: #FFFFFF;
          }
          .m3-gram-table td {
            min-width: 52px;
            text-align: center;
            padding: 5px 7px;
            border-radius: 5px;
            font-weight: 700;
            font-size: 0.86rem;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          }
          .m3-cell-ok { background: #DCFCE7; color: #166534; }
          .m3-cell-warn { background: #FEF3C7; color: #92400E; }
          .m3-cell-bad { background: #FEE2E2; color: #991B1B; }
          .m3-cell-neu { background: #F1F5F9; color: #334155; }
          .m3-badge {
            display: inline-block;
            padding: 5px 11px;
            border-radius: 999px;
            font-size: 0.82rem;
            font-weight: 700;
          }
          .m3-badge.ok { background: #DCFCE7; color: #166534; border: 1px solid #86EFAC; }
          .m3-badge.warn { background: #FEF3C7; color: #92400E; border: 1px solid #FCD34D; }
          .m3-badge.bad { background: #FEE2E2; color: #991B1B; border: 1px solid #FCA5A5; }
          .m3-plotly {
            width: 100%;
            height: 360px;
            border: 1px solid var(--border);
            border-radius: 10px;
            background: #FFFFFF;
            overflow: hidden;
          }
          .m3-iso-table {
            width: 100%;
            border-collapse: collapse;
            font-size: 0.86rem;
            margin-top: 10px;
          }
          .m3-iso-table th, .m3-iso-table td {
            border: 1px solid var(--border);
            padding: 7px 10px;
            text-align: center;
          }
          .m3-iso-table th { background: #F8FAFC; color: #334155; font-weight: 700; }
          .m3-quiz-item {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 14px 16px;
            margin-bottom: 12px;
          }
          .m3-opt-group { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 10px; }
          .m3-opt-btn {
            padding: 7px 13px;
            border: 1px solid #CBD5E1;
            border-radius: 8px;
            background: #FFFFFF;
            cursor: pointer;
            font-weight: 600;
            font-size: 0.88rem;
          }
          .m3-opt-btn:hover { background: #F1F5F9; }
          .m3-opt-btn.correct { background: #DCFCE7; border-color: #16A34A; color: #166534; }
          .m3-opt-btn.wrong { background: #FEE2E2; border-color: #DC2626; color: #991B1B; }
          .m3-fb { margin-top: 8px; font-size: 0.88rem; font-weight: 500; }
        </style>

        <!-- HEADER -->
        <div class="module-header">
          <span class="slide-ref">Part 3 · Lecture Notes ${meta.slides}</span>
          <h2>3. Orthogonal &amp; Orthonormal Sets, and Orthogonal Matrices</h2>
          <p>Why mutually perpendicular unit vectors make linear algebra <strong>effortless</strong>: automatic linear independence, instant coordinates without Gaussian elimination, and rigid transformations ($Q^T Q = I$).</p>
        </div>

        <!-- SECTION A: CORE DEFINITIONS & THEOREMS -->
        <div class="card">
          <h3>Orthogonal vs. Orthonormal Sets &amp; Effortless Coordinates</h3>
          <div class="grid-2">
            <div class="callout definition">
              <strong>Definitions:</strong><br/>
              A set of vectors $\\{\\mathbf{u}_1, \\dots, \\mathbf{u}_k\\}$ in $\\mathbb{R}^n$ is:
              <ul style="margin:6px 0 0 18px; padding:0; line-height:1.6;">
                <li><strong>Orthogonal</strong> if $\\langle \\mathbf{u}_i, \\mathbf{u}_j \\rangle = 0$ whenever $i \\ne j$.</li>
                <li><strong>Orthonormal</strong> if it is orthogonal <em>and</em> every vector has unit length $\\|\\mathbf{u}_i\\| = 1$, i.e. $\\langle \\mathbf{u}_i, \\mathbf{u}_j \\rangle = \\delta_{ij}$.</li>
              </ul>
            </div>

            <div class="callout tip">
              <strong>Theorem — Automatic Linear Independence:</strong><br/>
              Any orthogonal set of <strong>nonzero</strong> vectors $\\{\\mathbf{u}_1, \\dots, \\mathbf{u}_k\\}$ is linearly independent.<br/>
              <span style="font-size:0.86rem; color:var(--text-muted);"><strong>2-Line Proof:</strong> If $c_1\\mathbf{u}_1 + \\dots + c_k\\mathbf{u}_k = \\mathbf{0}$, take the inner product with $\\mathbf{u}_i$: all cross terms vanish, leaving $c_i\\|\\mathbf{u}_i\\|^2 = 0 \\implies c_i = 0$ (since $\\mathbf{u}_i \\ne \\mathbf{0}$).</span>
            </div>
          </div>

          <div class="callout definition" style="margin-top:12px;">
            <strong>Theorem — Coordinates Without Row Reduction!</strong><br/>
            If $\\{\\mathbf{u}_1, \\dots, \\mathbf{u}_k\\}$ is an <strong>orthogonal basis</strong> for a subspace $W$, then every $\\mathbf{x} \\in W$ expands uniquely as $\\mathbf{x} = c_1\\mathbf{u}_1 + \\dots + c_k\\mathbf{u}_k$, where:
            $$c_i = \\frac{\\langle \\mathbf{x}, \\mathbf{u}_i \\rangle}{\\langle \\mathbf{u}_i, \\mathbf{u}_i \\rangle} = \\frac{\\langle \\mathbf{x}, \\mathbf{u}_i \\rangle}{\\|\\mathbf{u}_i\\|^2}.$$
            If the basis $\\{\\mathbf{q}_1, \\dots, \\mathbf{q}_k\\}$ is <strong>orthonormal</strong>, the denominator is $1$, so simply $c_i = \\langle \\mathbf{x}, \\mathbf{q}_i \\rangle = \\mathbf{q}_i^T \\mathbf{x}$, and $\\|\\mathbf{x}\\|^2 = c_1^2 + \\dots + c_k^2$ (<em>Parseval's Identity</em>).
          </div>
        </div>

        <!-- SECTION B: INTERACTIVE LAB 1 — 2D ORTHONORMAL COORDINATES -->
        <div class="card">
          <h3>Interactive Lab 1: Coordinates in a Rotated Orthonormal Basis</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            Drag <strong style="color:#1D4ED8;">$\\mathbf{x}$ (blue)</strong> and spin the orthonormal basis $\\{\\mathbf{q}_1, \\mathbf{q}_2\\}$ with the angle slider $\\phi$. Notice how the coordinates $c_1, c_2$ are just dot products—and $c_1^2 + c_2^2 = \\|\\mathbf{x}\\|^2$ stays <strong>strictly invariant</strong>!
          </p>

          <div class="m3-controls">
            <label style="font-size:0.86rem; font-weight:700;">Basis Angle $\\phi$: <span id="m3PhiLbl" style="color:var(--polyu-red);">36.9°</span></label>
            <input type="range" id="m3PhiSlider" min="-180" max="180" step="0.1" value="36.87" style="width:170px;" />
            <button class="m3-chip" data-phi="0">Standard ($0^\\circ$)</button>
            <button class="m3-chip" data-phi="36.87">Rational $3\\text{-}4\\text{-}5$ ($36.9^\\circ$)</button>
            <button class="m3-chip" data-phi="45">$45^\\circ$ ($\\frac{1}{\\sqrt{2}}$)</button>
            <button class="m3-chip" data-phi="120">$120^\\circ$</button>
          </div>

          <div class="grid-2" style="align-items:start;">
            <div class="m3-canvas-wrap">
              <canvas id="m3CoordCanvas" class="m3-canvas"></canvas>
            </div>

            <div>
              <div class="m3-stat-grid">
                <div class="m3-stat">
                  <div class="m3-stat-label">Basis Vector $\\mathbf{q}_1$</div>
                  <div class="m3-stat-val" id="m3Q1Val" style="color:#0F766E;"></div>
                </div>
                <div class="m3-stat">
                  <div class="m3-stat-label">Basis Vector $\\mathbf{q}_2$</div>
                  <div class="m3-stat-val" id="m3Q2Val" style="color:#B45309;"></div>
                </div>
                <div class="m3-stat">
                  <div class="m3-stat-label">Coord $c_1 = \\langle \\mathbf{x}, \\mathbf{q}_1 \\rangle$</div>
                  <div class="m3-stat-val" id="m3C1Val" style="color:#0F766E;"></div>
                </div>
                <div class="m3-stat">
                  <div class="m3-stat-label">Coord $c_2 = \\langle \\mathbf{x}, \\mathbf{q}_2 \\rangle$</div>
                  <div class="m3-stat-val" id="m3C2Val" style="color:#B45309;"></div>
                </div>
              </div>

              <div class="callout tip" style="margin-top:12px;" id="m3CoordMath"></div>
            </div>
          </div>
        </div>

        <!-- SECTION C: ORTHOGONAL MATRICES & LAB 2 (GRAM MATRIX DIAGNOSTIC) -->
        <div class="card">
          <h3>Orthogonal Matrices ($Q^T Q = I$) &amp; Interactive Matrix Diagnostic</h3>
          <div class="grid-2" style="margin-bottom:14px;">
            <div class="callout definition">
              <strong>Matrices with Orthonormal Columns:</strong><br/>
              Let $U$ be an $m \\times n$ matrix with columns $\\mathbf{q}_1, \\dots, \\mathbf{q}_n$. Because entry $(i,j)$ of $U^T U$ is $\\mathbf{q}_i^T \\mathbf{q}_j = \\langle \\mathbf{q}_i, \\mathbf{q}_j \\rangle$:
              $$\\text{Columns of } U \\text{ are orthonormal} \\iff U^T U = I_n.$$
            </div>

            <div class="callout tip">
              <strong>Square Case — Orthogonal Matrix ($m = n$):</strong><br/>
              A <em>square</em> $n \\times n$ matrix $Q$ with orthonormal columns is called an <strong>orthogonal matrix</strong>:
              $$Q^T Q = I_n \\iff Q^{-1} = Q^T \\iff Q Q^T = I_n.$$
              <span style="font-size:0.84rem; color:var(--text-muted);"><strong>Watch out!</strong> For a tall $m \\times n$ matrix ($m > n$), $U^T U = I_n$ holds, but $U U^T \\ne I_m$! Test the $3\\times 2$ preset below.</span>
            </div>
          </div>

          <div class="m3-controls">
            <button class="m3-chip active" id="m3Mode33">Square $3 \\times 3$</button>
            <button class="m3-chip" id="m3Mode32">Tall $3 \\times 2$</button>
            <span style="color:#CBD5E1;">|</span>
            <button class="btn" id="m3PreRot">3D Rotation ($\\det=+1$)</button>
            <button class="btn" id="m3PreRef">Householder Reflection ($\\det=-1$)</button>
            <button class="btn" id="m3PreOrthoOnly">Orthogonal (Not Unit!)</button>
            <button class="btn" id="m3PreShear">Non-Orthogonal ($\\det=+1$)</button>
            <button class="btn" id="m3PreTall">Tall $3\\times 2$ Orthonormal</button>
          </div>

          <div class="grid-2" style="align-items:start;">
            <div>
              <div style="display:flex; align-items:center; gap:12px; flex-wrap:wrap; margin-bottom:12px;">
                <strong style="font-size:1.05rem;">$A =$</strong>
                <div id="m3MatInputWrap" class="m3-mat-grid"></div>
                <button class="m3-chip" id="m3NormalizeBtn" title="Divide each column by its Euclidean norm">Normalize Columns</button>
              </div>

              <div id="m3DiagBadge" style="margin-bottom:12px;"></div>

              <div style="display:flex; gap:18px; flex-wrap:wrap; align-items:flex-start;">
                <div>
                  <div class="m3-stat-label" style="margin-bottom:5px;">Gram Matrix $A^T A$</div>
                  <div id="m3GramATA"></div>
                </div>
                <div>
                  <div class="m3-stat-label" style="margin-bottom:5px;">Outer Product $A A^T$</div>
                  <div id="m3GramAAT"></div>
                </div>
              </div>

              <div id="m3DiagNotes" style="margin-top:10px; font-size:0.88rem; color:#334155;"></div>
            </div>

            <div id="m3Plotly3D" class="m3-plotly"></div>
          </div>

          <details class="code-details">
            <summary>Python (SymPy) Companion Code — Checking Orthogonality &amp; Computing $Q^{-1} = Q^T$</summary>
            <pre><code>import sympy as sp

# Exact rational 3x3 orthogonal matrix (Householder reflection)
Q = sp.Matrix([
    [1, 0, 0],
    [0, sp.Rational(7, 25), sp.Rational(-24, 25)],
    [0, sp.Rational(-24, 25), sp.Rational(-7, 25)]
])

print("Q^T * Q =", Q.T * Q)         # Identity matrix I_3
print("det(Q)  =", Q.det())         # -1 (orientation-reversing reflection)
print("Q^-1 == Q^T?", Q.inv() == Q.T)  # True!

# Tall 3x2 matrix with orthonormal columns
U = sp.Matrix([
    [sp.Rational(3, 5), sp.Rational(-12, 25)],
    [sp.Rational(4, 5), sp.Rational(9, 25)],
    [0,                 sp.Rational(4, 5)]
])
print("U^T * U =", U.T * U)         # I_2
print("U * U^T == I_3?", U * U.T == sp.eye(3))  # False! (It's a projection matrix!)</code></pre>
          </details>
        </div>

        <!-- SECTION D: INTERACTIVE LAB 3 — ISOMETRY / RIGID MOTION EXPLORER -->
        <div class="card">
          <h3>Interactive Lab 3: Isometries — Why Orthogonal Matrices Preserve Geometry</h3>
          <div class="callout definition" style="margin-bottom:12px;">
            <strong>Isometry Theorem:</strong> If $Q^T Q = I$ (columns of $Q$ are orthonormal), then for all $\\mathbf{u}, \\mathbf{v}$:
            $$\\langle Q\\mathbf{u}, Q\\mathbf{v} \\rangle = (Q\\mathbf{v})^T(Q\\mathbf{u}) = \\mathbf{v}^T (Q^T Q)\\mathbf{u} = \\mathbf{v}^T \\mathbf{u} = \\langle \\mathbf{u}, \\mathbf{v} \\rangle.$$
            Consequently, $Q$ preserves <strong>lengths</strong> ($\\|Q\\mathbf{u}\\| = \\|\\mathbf{u}\\|$), <strong>distances</strong> ($\\|Q\\mathbf{u} - Q\\mathbf{v}\\| = \\|\\mathbf{u} - \\mathbf{v}\\|$), and <strong>angles</strong>.
          </div>

          <div class="m3-controls">
            <button class="m3-chip active" data-iso="rot">Rotation $R_\\theta$ ($\\det = +1$)</button>
            <button class="m3-chip" data-iso="ref">Reflection $H_\\theta$ ($\\det = -1$)</button>
            <button class="m3-chip" data-iso="shear">Non-Orthogonal Shear ($\\det = +1$)</button>
            <label style="font-size:0.86rem; font-weight:700; margin-left:6px;">
              <span id="m3IsoParamName">Angle $\\theta$</span>: <span id="m3IsoParamVal" style="color:var(--polyu-red);">50°</span>
            </label>
            <input type="range" id="m3IsoSlider" min="-180" max="180" step="1" value="50" style="width:170px;" />
          </div>

          <div class="grid-2" style="align-items:start;">
            <div class="m3-canvas-wrap">
              <canvas id="m3IsoCanvas" class="m3-canvas-static"></canvas>
            </div>

            <div>
              <div id="m3IsoMatrixTex" style="background:#F8FAFC; border:1px solid var(--border); border-radius:8px; padding:10px 14px;"></div>
              <table class="m3-iso-table" id="m3IsoTable">
                <thead>
                  <tr>
                    <th>Quantity</th>
                    <th>Original ($\\mathbf{u}, \\mathbf{v}$)</th>
                    <th>Transformed ($M\\mathbf{u}, M\\mathbf{v}$)</th>
                    <th>Preserved?</th>
                  </tr>
                </thead>
                <tbody></tbody>
              </table>
              <p id="m3IsoComment" style="font-size:0.86rem; color:#334155; margin-top:10px;"></p>
            </div>
          </div>
        </div>

        <!-- SECTION E: QUICK CONCEPT CHECK -->
        <div class="card">
          <h3>Quick Concept Check (Exam Pitfall Radar)</h3>

          <div class="m3-quiz-item" data-q="1" data-ans="B">
            <strong>Q1.</strong> True or False: <em>"Every orthogonal set of vectors in $\\mathbb{R}^n$ is linearly independent."</em>
            <div class="m3-opt-group">
              <button class="m3-opt-btn" data-opt="A">True — orthogonality always guarantees linear independence.</button>
              <button class="m3-opt-btn" data-opt="B">False — it fails if the set contains the zero vector $\\mathbf{0}$.</button>
              <button class="m3-opt-btn" data-opt="C">False — it is only true for orthonormal sets.</button>
            </div>
            <div class="m3-fb" id="m3Fb1"></div>
          </div>

          <div class="m3-quiz-item" data-q="2" data-ans="A">
            <strong>Q2.</strong> Let $\\mathbf{u}_1 = \\begin{bmatrix}1\\\\1\\\\0\\end{bmatrix}$, $\\mathbf{u}_2 = \\begin{bmatrix}-1\\\\1\\\\0\\end{bmatrix}$, $\\mathbf{u}_3 = \\begin{bmatrix}0\\\\0\\\\2\\end{bmatrix}$ be an orthogonal basis of $\\mathbb{R}^3$. When writing $\\mathbf{x} = \\begin{bmatrix}3\\\\5\\\\4\\end{bmatrix} = c_1\\mathbf{u}_1 + c_2\\mathbf{u}_2 + c_3\\mathbf{u}_3$, what is $c_2$?
            <div class="m3-opt-group">
              <button class="m3-opt-btn" data-opt="A">$c_2 = 1$</button>
              <button class="m3-opt-btn" data-opt="B">$c_2 = 2$</button>
              <button class="m3-opt-btn" data-opt="C">$c_2 = \\sqrt{2}$</button>
              <button class="m3-opt-btn" data-opt="D">$c_2 = 4$</button>
            </div>
            <div class="m3-fb" id="m3Fb2"></div>
          </div>

          <div class="m3-quiz-item" data-q="3" data-ans="C">
            <strong>Q3.</strong> Suppose $U$ is a $4 \\times 2$ matrix with <strong>orthonormal columns</strong>. Which statement is ALWAYS true?
            <div class="m3-opt-group">
              <button class="m3-opt-btn" data-opt="A">$U U^T = I_4$</button>
              <button class="m3-opt-btn" data-opt="B">$U^{-1} = U^T$</button>
              <button class="m3-opt-btn" data-opt="C">$U^T U = I_2$</button>
              <button class="m3-opt-btn" data-opt="D">$\\det(U) = \\pm 1$</button>
            </div>
            <div class="m3-fb" id="m3Fb3"></div>
          </div>

          <div class="m3-quiz-item" data-q="4" data-ans="D">
            <strong>Q4.</strong> Let $Q$ be ANY $3 \\times 3$ orthogonal matrix, and let $\\mathbf{x} = \\begin{bmatrix} 3 \\\\ -4 \\\\ 0 \\end{bmatrix}$. What is $\\|Q\\mathbf{x}\\|$?
            <div class="m3-opt-group">
              <button class="m3-opt-btn" data-opt="A">Cannot be determined without knowing $Q$.</button>
              <button class="m3-opt-btn" data-opt="B">$25$</button>
              <button class="m3-opt-btn" data-opt="C">$\\pm 5$</button>
              <button class="m3-opt-btn" data-opt="D">$5$</button>
            </div>
            <div class="m3-fb" id="m3Fb4"></div>
          </div>
        </div>
      `;
    },

    mount(container, { renderMath }) {
      // =====================================================================
      // LAB 1: 2D ROTATED ORTHONORMAL COORDINATES CANVAS
      // =====================================================================
      const cCanvas = container.querySelector('#m3CoordCanvas');
      const cCtx = cCanvas.getContext('2d');
      const phiSlider = container.querySelector('#m3PhiSlider');
      const phiLbl = container.querySelector('#m3PhiLbl');

      let vx = [3.0, 2.0];
      let phiDeg = 36.87;
      let draggingX = false;
      const cRange = 5.2;

      function cLayout() {
        const r = cCanvas.getBoundingClientRect();
        const w = r.width || 450, h = r.height || 360;
        const sc = Math.min(w, h) / (2 * cRange);
        return { w, h, cx: w / 2, cy: h / 2, sc };
      }
      function cToScr(pt, L) { return [L.cx + pt[0] * L.sc, L.cy - pt[1] * L.sc]; }
      function cToWld(sx, sy, L) {
        let wx = Math.round(((sx - L.cx) / L.sc) * 4) / 4;
        let wy = Math.round(((L.cy - sy) / L.sc) * 4) / 4;
        return [Math.max(-4.5, Math.min(4.5, wx)), Math.max(-4.5, Math.min(4.5, wy))];
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

      function drawLab1() {
        if (!cCanvas.isConnected) return;
        const dpr = window.devicePixelRatio || 1;
        const L = cLayout();
        cCanvas.width = Math.round(L.w * dpr);
        cCanvas.height = Math.round(L.h * dpr);
        cCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
        cCtx.clearRect(0, 0, L.w, L.h);

        const rad = (phiDeg * Math.PI) / 180;
        const q1 = [Math.cos(rad), Math.sin(rad)];
        const q2 = [-Math.sin(rad), Math.cos(rad)];

        // Faint standard axes
        cCtx.strokeStyle = '#E2E8F0';
        cCtx.lineWidth = 1.2;
        cCtx.beginPath(); cCtx.moveTo(0, L.cy); cCtx.lineTo(L.w, L.cy); cCtx.stroke();
        cCtx.beginPath(); cCtx.moveTo(L.cx, 0); cCtx.lineTo(L.cx, L.h); cCtx.stroke();

        // Rotated orthonormal grid lines
        cCtx.save();
        cCtx.lineWidth = 1;
        for (let k = -7; k <= 7; k++) {
          const isAxis = k === 0;
          cCtx.strokeStyle = isAxis ? 'rgba(15, 118, 110, 0.35)' : 'rgba(15, 118, 110, 0.11)';
          const pA = cToScr(add(scale(-9, q1), scale(k, q2)), L);
          const pB = cToScr(add(scale(9, q1), scale(k, q2)), L);
          cCtx.beginPath(); cCtx.moveTo(pA[0], pA[1]); cCtx.lineTo(pB[0], pB[1]); cCtx.stroke();

          cCtx.strokeStyle = isAxis ? 'rgba(180, 83, 9, 0.35)' : 'rgba(180, 83, 9, 0.11)';
          const rA = cToScr(add(scale(k, q1), scale(-9, q2)), L);
          const rB = cToScr(add(scale(k, q1), scale(9, q2)), L);
          cCtx.beginPath(); cCtx.moveTo(rA[0], rA[1]); cCtx.lineTo(rB[0], rB[1]); cCtx.stroke();
        }
        cCtx.restore();

        const c1 = dot(vx, q1);
        const c2 = dot(vx, q2);
        const comp1 = scale(c1, q1);
        const comp2 = scale(c2, q2);

        const [xS, yS] = cToScr(vx, L);
        const [p1x, p1y] = cToScr(comp1, L);
        const [p2x, p2y] = cToScr(comp2, L);

        // Dashed rectangle showing orthogonal decomposition
        cCtx.save();
        cCtx.strokeStyle = '#94A3B8';
        cCtx.lineWidth = 1.6;
        cCtx.setLineDash([4, 4]);
        cCtx.beginPath(); cCtx.moveTo(p1x, p1y); cCtx.lineTo(xS, yS); cCtx.stroke();
        cCtx.beginPath(); cCtx.moveTo(p2x, p2y); cCtx.lineTo(xS, yS); cCtx.stroke();
        cCtx.restore();

        // Components c1*q1 and c2*q2
        if (Math.abs(c1) > 0.05) drawArrow(cCtx, L.cx, L.cy, p1x, p1y, 'rgba(15,118,110,0.45)', 4.5);
        if (Math.abs(c2) > 0.05) drawArrow(cCtx, L.cx, L.cy, p2x, p2y, 'rgba(217,119,6,0.45)', 4.5);

        // Unit circle
        cCtx.save();
        cCtx.strokeStyle = 'rgba(100, 116, 139, 0.3)';
        cCtx.setLineDash([3, 3]);
        cCtx.beginPath();
        cCtx.arc(L.cx, L.cy, L.sc, 0, Math.PI * 2);
        cCtx.stroke();
        cCtx.restore();

        // Unit basis vectors q1, q2
        const [q1x, q1y] = cToScr(q1, L);
        const [q2x, q2y] = cToScr(q2, L);
        drawArrow(cCtx, L.cx, L.cy, q1x, q1y, '#0F766E', 3.0);
        drawArrow(cCtx, L.cx, L.cy, q2x, q2y, '#B45309', 3.0);

        cCtx.font = 'bold 13px sans-serif';
        cCtx.fillStyle = '#0F766E';
        cCtx.fillText('q₁', q1x + 7, q1y - 6);
        cCtx.fillStyle = '#B45309';
        cCtx.fillText('q₂', q2x + 7, q2y - 6);

        // Vector x
        drawArrow(cCtx, L.cx, L.cy, xS, yS, '#1D4ED8', 3.2);
        cCtx.beginPath();
        cCtx.arc(xS, yS, 7, 0, Math.PI * 2);
        cCtx.fillStyle = '#FFFFFF';
        cCtx.fill();
        cCtx.lineWidth = 2.5;
        cCtx.strokeStyle = '#1D4ED8';
        cCtx.stroke();
        cCtx.fillStyle = '#1D4ED8';
        cCtx.font = 'bold 14px sans-serif';
        cCtx.fillText('x', xS + 10, yS - 8);

        // Update readouts
        phiLbl.textContent = `${phiDeg.toFixed(1)}°`;
        container.querySelector('#m3Q1Val').textContent = `(${q1[0].toFixed(2)}, ${q1[1].toFixed(2)})ᵀ`;
        container.querySelector('#m3Q2Val').textContent = `(${q2[0].toFixed(2)}, ${q2[1].toFixed(2)})ᵀ`;
        container.querySelector('#m3C1Val').textContent = c1.toFixed(3);
        container.querySelector('#m3C2Val').textContent = c2.toFixed(3);

        const normSq = dot(vx, vx);
        const parseval = c1 * c1 + c2 * c2;
        const coordEl = container.querySelector('#m3CoordMath');
        coordEl.innerHTML = `
          <strong>Live Reconstruction &amp; Parseval's Identity:</strong><br/>
          $\\mathbf{x} = \\begin{bmatrix} ${vx[0].toFixed(2)} \\\\ ${vx[1].toFixed(2)} \\end{bmatrix} = (${c1.toFixed(2)})\\mathbf{q}_1 + (${c2.toFixed(2)})\\mathbf{q}_2$<br/>
          <span style="color:#15803D; font-weight:600;">
            $c_1^2 + c_2^2 = ${(c1 * c1).toFixed(2)} + ${(c2 * c2).toFixed(2)} = ${parseval.toFixed(2)} = \\|\\mathbf{x}\\|^2 = ${normSq.toFixed(2)}$ ✓
          </span>
        `;
        renderMath(coordEl);
      }

      function cPtr(e) {
        const r = cCanvas.getBoundingClientRect();
        const cx = e.touches ? e.touches[0].clientX : e.clientX;
        const cy = e.touches ? e.touches[0].clientY : e.clientY;
        return [cx - r.left, cy - r.top];
      }
      cCanvas.addEventListener('mousedown', e => {
        const L = cLayout(), [px, py] = cPtr(e), [xs, ys] = cToScr(vx, L);
        if (Math.hypot(px - xs, py - ys) < 24) { draggingX = true; e.preventDefault(); }
      });
      cCanvas.addEventListener('touchstart', e => {
        const L = cLayout(), [px, py] = cPtr(e), [xs, ys] = cToScr(vx, L);
        if (Math.hypot(px - xs, py - ys) < 28) { draggingX = true; e.preventDefault(); }
      }, { passive: false });
      window.addEventListener('mousemove', e => {
        if (!draggingX || !cCanvas.isConnected) return;
        const L = cLayout(), [px, py] = cPtr(e);
        vx = cToWld(px, py, L);
        drawLab1();
      });
      window.addEventListener('touchmove', e => {
        if (!draggingX || !cCanvas.isConnected) return;
        e.preventDefault();
        const L = cLayout(), [px, py] = cPtr(e);
        vx = cToWld(px, py, L);
        drawLab1();
      }, { passive: false });
      window.addEventListener('mouseup', () => { draggingX = false; });
      window.addEventListener('touchend', () => { draggingX = false; });

      phiSlider.addEventListener('input', () => {
        phiDeg = parseFloat(phiSlider.value);
        drawLab1();
      });
      container.querySelectorAll('[data-phi]').forEach(btn => {
        btn.addEventListener('click', () => {
          phiDeg = parseFloat(btn.dataset.phi);
          phiSlider.value = phiDeg;
          drawLab1();
        });
      });

      // =====================================================================
      // LAB 2: GRAM MATRIX DIAGNOSTIC (3x3 vs 3x2) + 3D PLOTLY FRAME
      // =====================================================================
      let matCols = 3;
      let matData = [
        [0.6, -0.8, 0],
        [0.8, 0.6, 0],
        [0, 0, 1]
      ];

      const matWrap = container.querySelector('#m3MatInputWrap');
      const mode33Btn = container.querySelector('#m3Mode33');
      const mode32Btn = container.querySelector('#m3Mode32');
      const plotly3D = container.querySelector('#m3Plotly3D');

      function buildMatrixInputs() {
        matWrap.style.gridTemplateColumns = `repeat(${matCols}, 68px)`;
        matWrap.innerHTML = '';
        for (let r = 0; r < 3; r++) {
          for (let c = 0; c < matCols; c++) {
            const inp = document.createElement('input');
            inp.type = 'number';
            inp.step = 'any';
            inp.className = 'm3-mat-input';
            inp.value = cleanNum(matData[r][c], 4);
            inp.dataset.r = r;
            inp.dataset.c = c;
            inp.addEventListener('input', () => {
              const v = parseFloat(inp.value);
              matData[r][c] = Number.isFinite(v) ? v : 0;
              updateLab2();
            });
            matWrap.appendChild(inp);
          }
        }
      }

      function renderGramTable(M, checkIdentity) {
        const rows = M.length, cols = M[0].length;
        let html = '<table class="m3-gram-table">';
        for (let i = 0; i < rows; i++) {
          html += '<tr>';
          for (let j = 0; j < cols; j++) {
            const v = M[i][j];
            let cls = 'm3-cell-neu';
            if (checkIdentity) {
              if (i === j) {
                if (Math.abs(v - 1) < 1e-3) cls = 'm3-cell-ok';
                else if (v > 1e-6) cls = 'm3-cell-warn';
                else cls = 'm3-cell-bad';
              } else {
                cls = Math.abs(v) < 1e-3 ? 'm3-cell-ok' : 'm3-cell-bad';
              }
            } else {
              const isI = (i === j && Math.abs(v - 1) < 1e-3) || (i !== j && Math.abs(v) < 1e-3);
              cls = isI ? 'm3-cell-ok' : 'm3-cell-neu';
            }
            html += `<td class="${cls}">${cleanNum(v, 2)}</td>`;
          }
          html += '</tr>';
        }
        html += '</table>';
        return html;
      }

      function updateLab2() {
        const A = matData.map(row => row.slice(0, matCols));
        const AT = transpose(A);
        const ATA = matMul(AT, A);
        const AAT = matMul(A, AT);

        // Classify columns
        let allOffZero = true, allDiagOne = true, hasZeroCol = false;
        for (let i = 0; i < matCols; i++) {
          if (ATA[i][i] < 1e-6) hasZeroCol = true;
          if (Math.abs(ATA[i][i] - 1) > 1e-3) allDiagOne = false;
          for (let j = 0; j < matCols; j++) {
            if (i !== j && Math.abs(ATA[i][j]) > 1e-3) allOffZero = false;
          }
        }

        const badgeEl = container.querySelector('#m3DiagBadge');
        const notesEl = container.querySelector('#m3DiagNotes');
        let noteHtml = '';

        if (allOffZero && allDiagOne) {
          if (matCols === 3) {
            const d = det3(A);
            const motion = d > 0 ? 'Proper Rotation (preserves orientation)' : 'Reflection / Improper Rotation (reverses orientation)';
            badgeEl.innerHTML = `<span class="m3-badge ok">✓ Orthogonal Matrix ($Q^T Q = Q Q^T = I_3$)</span>`;
            noteHtml = `All 3 columns are orthonormal! $\\det(A) = \\mathbf{${cleanNum(d, 3)}}$ (${motion}). Inverse is simply $A^{-1} = A^T$.`;
          } else {
            badgeEl.innerHTML = `<span class="m3-badge ok">✓ Orthonormal Columns ($U^T U = I_2$)</span>`;
            noteHtml = `Notice: $A^T A = I_2$, <strong>but $A A^T \\ne I_3$</strong>! Instead, $A A^T$ is the orthogonal projection matrix onto $\\operatorname{Col}(A)$ (preview of Module 4!).`;
          }
        } else if (allOffZero && !hasZeroCol) {
          badgeEl.innerHTML = `<span class="m3-badge warn">⚠ Orthogonal Columns (Not Unit Length)</span>`;
          noteHtml = `Off-diagonals of $A^T A$ are $0$, so columns are mutually perpendicular—but norms $\\ne 1$. Click <strong>Normalize Columns</strong> to make $A^T A = I$!`;
        } else {
          badgeEl.innerHTML = `<span class="m3-badge bad">✗ Not Orthogonal (Off-Diagonal $\\ne 0$)</span>`;
          const detStr = matCols === 3 ? ` (Even though $\\det(A) = ${cleanNum(det3(A), 2)}$, non-zero off-diagonals in $A^T A$ prove the columns are not perpendicular!)` : '';
          noteHtml = `At least one pair of columns has $\\langle \\mathbf{a}_i, \\mathbf{a}_j \\rangle \\ne 0$.${detStr}`;
        }

        container.querySelector('#m3GramATA').innerHTML = renderGramTable(ATA, true);
        container.querySelector('#m3GramAAT').innerHTML = renderGramTable(AAT, false);
        notesEl.innerHTML = noteHtml;
        renderMath(badgeEl);
        renderMath(notesEl);

        // 3D Plotly column visualization
        if (typeof Plotly === 'undefined') {
          plotly3D.innerHTML = '<div style="padding:20px; color:var(--text-muted);">3D plot requires Plotly.js.</div>';
          return;
        }
        const colors = ['#1D4ED8', '#DC2626', '#0F766E'];
        const traces = [];
        for (let c = 0; c < matCols; c++) {
          const v = [A[0][c], A[1][c], A[2][c]];
          traces.push({
            type: 'scatter3d',
            mode: 'lines+markers+text',
            name: `col ${c + 1} (‖a${c + 1}‖=${norm(v).toFixed(2)})`,
            x: [0, v[0]], y: [0, v[1]], z: [0, v[2]],
            text: ['', `a${c + 1}`],
            textposition: 'top center',
            textfont: { color: colors[c], size: 13 },
            line: { color: colors[c], width: 6 },
            marker: { size: [3, 6], color: colors[c] }
          });
        }
        // Origin reference point
        traces.push({
          type: 'scatter3d', mode: 'markers', showlegend: false,
          x: [0], y: [0], z: [0], marker: { size: 4, color: '#334155' }
        });

        Plotly.react(plotly3D, traces, {
          margin: { l: 0, r: 0, t: 10, b: 0 },
          legend: { orientation: 'h', y: 0.99, x: 0.02, font: { size: 11 } },
          scene: {
            aspectmode: 'data',
            xaxis: { title: 'x₁', zerolinecolor: '#94A3B8' },
            yaxis: { title: 'x₂', zerolinecolor: '#94A3B8' },
            zaxis: { title: 'x₃', zerolinecolor: '#94A3B8' },
            camera: { eye: { x: 1.5, y: 1.45, z: 1.1 } }
          }
        }, { responsive: true, displayModeBar: false });
      }

      function loadMatrixPreset(cols, data) {
        matCols = cols;
        mode33Btn.classList.toggle('active', cols === 3);
        mode32Btn.classList.toggle('active', cols === 2);
        matData = [0, 1, 2].map(r => [data[r][0], data[r][1], data[r][2] ?? 0]);
        buildMatrixInputs();
        updateLab2();
      }

      mode33Btn.addEventListener('click', () => loadMatrixPreset(3, matData));
      mode32Btn.addEventListener('click', () => loadMatrixPreset(2, matData));

      container.querySelector('#m3PreRot').addEventListener('click', () =>
        loadMatrixPreset(3, [[0.6, -0.8, 0], [0.8, 0.6, 0], [0, 0, 1]]));
      container.querySelector('#m3PreRef').addEventListener('click', () =>
        loadMatrixPreset(3, [[1, 0, 0], [0, 0.28, -0.96], [0, -0.96, -0.28]]));
      container.querySelector('#m3PreOrthoOnly').addEventListener('click', () =>
        loadMatrixPreset(3, [[1, 1, 1], [1, -1, 1], [1, 0, -2]]));
      container.querySelector('#m3PreShear').addEventListener('click', () =>
        loadMatrixPreset(3, [[1, 1, 0], [0, 1, 1], [0, 0, 1]]));
      container.querySelector('#m3PreTall').addEventListener('click', () =>
        loadMatrixPreset(2, [[0.6, -0.48], [0.8, 0.36], [0, 0.8]]));

      container.querySelector('#m3NormalizeBtn').addEventListener('click', () => {
        for (let c = 0; c < matCols; c++) {
          const n = norm([matData[0][c], matData[1][c], matData[2][c]]);
          if (n > 1e-8) {
            matData[0][c] = Number((matData[0][c] / n).toFixed(4));
            matData[1][c] = Number((matData[1][c] / n).toFixed(4));
            matData[2][c] = Number((matData[2][c] / n).toFixed(4));
          }
        }
        buildMatrixInputs();
        updateLab2();
      });

      // =====================================================================
      // LAB 3: ISOMETRY / RIGID MOTION VISUALIZER
      // =====================================================================
      const isoCanvas = container.querySelector('#m3IsoCanvas');
      const isoCtx = isoCanvas.getContext('2d');
      const isoSlider = container.querySelector('#m3IsoSlider');
      const isoParamName = container.querySelector('#m3IsoParamName');
      const isoParamVal = container.querySelector('#m3IsoParamVal');
      let isoMode = 'rot';

      // Two test vectors + asymmetric L-flag polygon to reveal orientation flips
      const u0 = [2.2, 0.6];
      const v0 = [0.7, 1.9];
      const flagPoly = [[0, 0], [2.2, 0.6], [1.7, 1.7], [0.7, 1.9]];

      function getIsoMatrix() {
        const val = parseFloat(isoSlider.value);
        if (isoMode === 'rot') {
          const r = (val * Math.PI) / 180;
          return {
            M: [[Math.cos(r), -Math.sin(r)], [Math.sin(r), Math.cos(r)]],
            det: 1,
            label: `${val}°`,
            tex: `Q = R_{${val}^\\circ} = \\begin{bmatrix} ${Math.cos(r).toFixed(2)} & ${(-Math.sin(r)).toFixed(2)} \\\\ ${Math.sin(r).toFixed(2)} & ${Math.cos(r).toFixed(2)} \\end{bmatrix}, \\quad \\det(Q) = +1`,
            comment: 'Rotation is an orthogonal matrix (QᵀQ = I, det = +1): every length, angle, and orientation is strictly preserved.'
          };
        } else if (isoMode === 'ref') {
          const r = (val * Math.PI) / 180; // reflection across line at angle r/2
          return {
            M: [[Math.cos(r), Math.sin(r)], [Math.sin(r), -Math.cos(r)]],
            det: -1,
            label: `${(val / 2).toFixed(1)}°`,
            tex: `Q = H = \\begin{bmatrix} ${Math.cos(r).toFixed(2)} & ${Math.sin(r).toFixed(2)} \\\\ ${Math.sin(r).toFixed(2)} & ${(-Math.cos(r)).toFixed(2)} \\end{bmatrix}, \\quad \\det(Q) = -1`,
            comment: 'Reflection is an orthogonal matrix (QᵀQ = I, det = -1): lengths and angles are preserved, while orientation flips (mirror image).'
          };
        } else {
          const k = val / 90; // shear factor in [-2, 2]
          return {
            M: [[1, k], [0, 1]],
            det: 1,
            label: `k = ${k.toFixed(2)}`,
            tex: `M = \\begin{bmatrix} 1 & ${k.toFixed(2)} \\\\ 0 & 1 \\end{bmatrix}, \\quad \\det(M) = +1 \\text{ (NOT orthogonal!)}`,
            comment: 'Notice! Even though det(M) = +1 (area is preserved), MᵀM ≠ I, so lengths, inner products, and angles are distorted!'
          };
        }
      }

      function drawLab3() {
        if (!isoCanvas.isConnected) return;
        const { M, label, tex, comment } = getIsoMatrix();
        isoParamVal.textContent = label;

        const dpr = window.devicePixelRatio || 1;
        const rect = isoCanvas.getBoundingClientRect();
        const w = rect.width || 440, h = rect.height || 320;
        isoCanvas.width = Math.round(w * dpr);
        isoCanvas.height = Math.round(h * dpr);
        isoCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
        isoCtx.clearRect(0, 0, w, h);

        const cx = w / 2, cy = h / 2, sc = Math.min(w, h) / 9.2;
        const toS = p => [cx + p[0] * sc, cy - p[1] * sc];

        // Grid + axes
        isoCtx.strokeStyle = '#F1F5F9'; isoCtx.lineWidth = 1;
        for (let g = -5; g <= 5; g++) {
          const [gx] = toS([g, 0]), [, gy] = toS([0, g]);
          isoCtx.beginPath(); isoCtx.moveTo(gx, 0); isoCtx.lineTo(gx, h); isoCtx.stroke();
          isoCtx.beginPath(); isoCtx.moveTo(0, gy); isoCtx.lineTo(w, gy); isoCtx.stroke();
        }
        isoCtx.strokeStyle = '#CBD5E1'; isoCtx.lineWidth = 1.3;
        isoCtx.beginPath(); isoCtx.moveTo(0, cy); isoCtx.lineTo(w, cy); isoCtx.stroke();
        isoCtx.beginPath(); isoCtx.moveTo(cx, 0); isoCtx.lineTo(cx, h); isoCtx.stroke();

        // If reflection, draw mirror axis
        if (isoMode === 'ref') {
          const half = (parseFloat(isoSlider.value) * Math.PI) / 360;
          const a = toS([-6 * Math.cos(half), -6 * Math.sin(half)]);
          const b = toS([6 * Math.cos(half), 6 * Math.sin(half)]);
          isoCtx.save();
          isoCtx.strokeStyle = 'rgba(160, 35, 55, 0.38)';
          isoCtx.lineWidth = 2;
          isoCtx.setLineDash([6, 4]);
          isoCtx.beginPath(); isoCtx.moveTo(a[0], a[1]); isoCtx.lineTo(b[0], b[1]); isoCtx.stroke();
          isoCtx.restore();
        }

        // Draw Original polygon (dashed ghost)
        isoCtx.save();
        isoCtx.fillStyle = 'rgba(148, 163, 184, 0.16)';
        isoCtx.strokeStyle = '#94A3B8';
        isoCtx.lineWidth = 1.6;
        isoCtx.setLineDash([4, 4]);
        isoCtx.beginPath();
        flagPoly.forEach((p, idx) => {
          const [sx, sy] = toS(p);
          if (idx === 0) isoCtx.moveTo(sx, sy); else isoCtx.lineTo(sx, sy);
        });
        isoCtx.closePath(); isoCtx.fill(); isoCtx.stroke();
        isoCtx.restore();

        // Draw Transformed polygon
        const transPoly = flagPoly.map(p => matVec2(M, p));
        isoCtx.save();
        isoCtx.fillStyle = isoMode === 'shear' ? 'rgba(220, 38, 38, 0.15)' : 'rgba(15, 118, 110, 0.16)';
        isoCtx.strokeStyle = isoMode === 'shear' ? '#DC2626' : '#0F766E';
        isoCtx.lineWidth = 2.2;
        isoCtx.beginPath();
        transPoly.forEach((p, idx) => {
          const [sx, sy] = toS(p);
          if (idx === 0) isoCtx.moveTo(sx, sy); else isoCtx.lineTo(sx, sy);
        });
        isoCtx.closePath(); isoCtx.fill(); isoCtx.stroke();
        isoCtx.restore();

        // Original vectors (faint) & Transformed vectors (bold)
        const u1 = matVec2(M, u0), v1 = matVec2(M, v0);
        const [u0x, u0y] = toS(u0), [v0x, v0y] = toS(v0);
        const [u1x, u1y] = toS(u1), [v1x, v1y] = toS(v1);

        drawArrow(isoCtx, cx, cy, u0x, u0y, 'rgba(29,78,216,0.35)', 2.2, true);
        drawArrow(isoCtx, cx, cy, v0x, v0y, 'rgba(217,119,6,0.35)', 2.2, true);
        drawArrow(isoCtx, cx, cy, u1x, u1y, '#1D4ED8', 3.0);
        drawArrow(isoCtx, cx, cy, v1x, v1y, '#D97706', 3.0);

        isoCtx.font = 'bold 12px sans-serif';
        isoCtx.fillStyle = '#1D4ED8'; isoCtx.fillText('Mu', u1x + 6, u1y - 5);
        isoCtx.fillStyle = '#B45309'; isoCtx.fillText('Mv', v1x + 6, v1y - 5);

        // Metrics comparison table
        const nu0 = norm(u0), nu1 = norm(u1);
        const nv0 = norm(v0), nv1 = norm(v1);
        const d0 = dot(u0, v0), d1 = dot(u1, v1);
        const ang0 = (Math.acos(Math.max(-1, Math.min(1, d0 / (nu0 * nv0)))) * 180) / Math.PI;
        const ang1 = (Math.acos(Math.max(-1, Math.min(1, d1 / (nu1 * nv1)))) * 180) / Math.PI;

        const rows = [
          ['Length ‖u‖', nu0.toFixed(3), nu1.toFixed(3), Math.abs(nu0 - nu1) < 1e-3],
          ['Length ‖v‖', nv0.toFixed(3), nv1.toFixed(3), Math.abs(nv0 - nv1) < 1e-3],
          ['Dot ⟨u, v⟩', d0.toFixed(3), d1.toFixed(3), Math.abs(d0 - d1) < 1e-3],
          ['Angle θ', `${ang0.toFixed(1)}°`, `${ang1.toFixed(1)}°`, Math.abs(ang0 - ang1) < 0.1]
        ];
        const tbody = container.querySelector('#m3IsoTable tbody');
        tbody.innerHTML = rows.map(([name, a, b, ok]) => `
          <tr>
            <td><strong>${name}</strong></td>
            <td>${a}</td>
            <td>${b}</td>
            <td style="font-weight:700; color:${ok ? '#15803D' : '#DC2626'};">${ok ? '✓ Equal' : '✗ Changed'}</td>
          </tr>
        `).join('');

        const matEl = container.querySelector('#m3IsoMatrixTex');
        matEl.innerHTML = `$$${tex}$$`;
        renderMath(matEl);
        container.querySelector('#m3IsoComment').textContent = comment;
      }

      container.querySelectorAll('[data-iso]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-iso]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          isoMode = btn.dataset.iso;
          isoParamName.textContent = isoMode === 'rot' ? 'Rotation Angle' : isoMode === 'ref' ? 'Mirror Line Angle' : 'Shear Factor';
          drawLab3();
        });
      });
      isoSlider.addEventListener('input', drawLab3);

      // =====================================================================
      // SECTION E: QUICK CONCEPT CHECK
      // =====================================================================
      const explanations = {
        1: 'The zero vector 0 is orthogonal to every vector, yet any set containing 0 is automatically linearly dependent! The theorem requires NONZERO orthogonal vectors.',
        2: 'Since {u₁, u₂, u₃} is orthogonal (not orthonormal!), we MUST divide by ‖u₂‖² = (-1)² + 1² = 2. Thus c₂ = ⟨x, u₂⟩ / ‖u₂‖² = (-3 + 5) / 2 = 1.',
        3: 'For a 4×2 matrix U with orthonormal columns, entry (i,j) of the 2×2 matrix UᵀU is ⟨qᵢ, qⱼ⟩ = δᵢⱼ, so UᵀU = I₂. Meanwhile, UUᵀ is 4×4 with rank ≤ 2, so UUᵀ ≠ I₄.',
        4: 'Orthogonal matrices are isometries: ‖Qx‖ = ‖x‖ = √(3² + (-4)² + 0²) = 5 (norms are always non-negative!).'
      };

      container.querySelectorAll('.m3-quiz-item').forEach(item => {
        const qNum = item.dataset.q;
        const ans = item.dataset.ans;
        const fb = container.querySelector(`#m3Fb${qNum}`);
        item.querySelectorAll('.m3-opt-btn').forEach(btn => {
          btn.addEventListener('click', () => {
            item.querySelectorAll('.m3-opt-btn').forEach(b => b.classList.remove('correct', 'wrong'));
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

      window.addEventListener('resize', () => {
        drawLab1();
        drawLab3();
      });

      // Initial renders
      drawLab1();
      buildMatrixInputs();
      updateLab2();
      drawLab3();
    }
  });
})();
