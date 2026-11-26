/**
 * AMA1751 Linear Algebra — Topic 6: Inner Product
 * Module 4: Orthogonal Complements, Subspace Projections & Best Approximation
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
  function cross3(a, b) {
    return [
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0]
    ];
  }
  function matVec(M, v) {
    return M.map(row => dot(row, v));
  }
  function outerProd(u, v) {
    return u.map(ui => v.map(vj => ui * vj));
  }
  function addMat(A, B) {
    return A.map((row, i) => row.map((v, j) => v + B[i][j]));
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
  function clean(x, d = 2) {
    if (Math.abs(x) < 1e-6) return '0';
    return String(Number(x.toFixed(d)));
  }
  function fmtVec(v, d = 2) {
    return `(${v.map(x => clean(x, d)).join(', ')})ᵀ`;
  }
  // 2-vector or general Gram-Schmidt returning orthonormal basis
  function gramSchmidt(vecs) {
    const basis = [];
    for (const v of vecs) {
      let w = v.slice();
      for (const q of basis) {
        w = sub(w, scale(dot(w, q), q));
      }
      const n = norm(w);
      if (n > 1e-8) basis.push(scale(1 / n, w));
    }
    return basis;
  }
  function renderMatHTML(M, d = 2) {
    let html = '<table class="m4-mat-table">';
    for (let i = 0; i < M.length; i++) {
      html += '<tr>';
      for (let j = 0; j < M[0].length; j++) {
        html += `<td>${clean(M[i][j], d)}</td>`;
      }
      html += '</tr>';
    }
    return html + '</table>';
  }

  AMA1751.registerModule(4, {
    render(meta) {
      return `
        <style>
          .m4-canvas-wrap {
            position: relative;
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 10px;
            overflow: hidden;
          }
          .m4-canvas {
            display: block;
            width: 100%;
            height: 370px;
            cursor: grab;
            touch-action: none;
          }
          .m4-canvas:active { cursor: grabbing; }
          .m4-controls {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            margin-bottom: 12px;
            align-items: center;
          }
          .m4-chip {
            padding: 6px 11px;
            border-radius: 999px;
            border: 1px solid var(--border);
            background: #F8FAFC;
            font-size: 0.82rem;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.15s;
          }
          .m4-chip:hover, .m4-chip.active {
            background: var(--polyu-red-soft);
            border-color: #E8B4C4;
            color: var(--polyu-red-dark);
          }
          .m4-stat-grid {
            display: grid;
            grid-template-columns: repeat(2, minmax(0, 1fr));
            gap: 9px;
          }
          .m4-stat {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 9px 11px;
          }
          .m4-stat-label {
            font-size: 0.73rem;
            color: var(--text-muted);
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.04em;
          }
          .m4-stat-val {
            font-size: 0.94rem;
            font-weight: 700;
            margin-top: 3px;
            color: #0F172A;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          }
          .m4-mat-table {
            border-collapse: separate;
            border-spacing: 3px;
            display: inline-table;
            padding: 5px 9px;
            border-left: 2px solid #475569;
            border-right: 2px solid #475569;
            border-radius: 6px;
            background: #FFFFFF;
          }
          .m4-mat-table td {
            min-width: 46px;
            text-align: center;
            padding: 4px 5px;
            border-radius: 4px;
            background: #F8FAFC;
            font-weight: 700;
            font-size: 0.84rem;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          }
          .m4-plotly {
            width: 100%;
            height: 410px;
            border: 1px solid var(--border);
            border-radius: 10px;
            background: #FFFFFF;
            overflow: hidden;
          }
          .m4-slider-row {
            display: grid;
            grid-template-columns: 78px 1fr 44px;
            align-items: center;
            gap: 8px;
            margin-bottom: 6px;
            font-size: 0.85rem;
            font-weight: 600;
          }
          .m4-fs-grid {
            display: grid;
            grid-template-columns: repeat(2, minmax(0, 1fr));
            gap: 12px;
            margin-top: 12px;
          }
          @media (max-width: 720px) {
            .m4-fs-grid { grid-template-columns: 1fr; }
          }
          .m4-fs-box {
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 12px 14px;
            background: #F8FAFC;
          }
          .m4-fs-box h4 {
            margin: 0 0 6px 0;
            font-size: 0.92rem;
            display: flex;
            justify-content: space-between;
            align-items: center;
          }
          .m4-dim-pill {
            font-size: 0.75rem;
            padding: 2px 8px;
            border-radius: 999px;
            font-weight: 700;
          }
          .m4-num-input {
            width: 62px;
            padding: 5px 6px;
            border: 1px solid #CBD5E1;
            border-radius: 6px;
            text-align: center;
            font-weight: 700;
            font-size: 0.9rem;
          }
          .m4-quiz-item {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 14px 16px;
            margin-bottom: 12px;
          }
          .m4-opt-group { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 10px; }
          .m4-opt-btn {
            padding: 7px 13px;
            border: 1px solid #CBD5E1;
            border-radius: 8px;
            background: #FFFFFF;
            cursor: pointer;
            font-weight: 600;
            font-size: 0.88rem;
          }
          .m4-opt-btn:hover { background: #F1F5F9; }
          .m4-opt-btn.correct { background: #DCFCE7; border-color: #16A34A; color: #166534; }
          .m4-opt-btn.wrong { background: #FEE2E2; border-color: #DC2626; color: #991B1B; }
          .m4-fb { margin-top: 8px; font-size: 0.88rem; font-weight: 500; }
        </style>

        <!-- HEADER -->
        <div class="module-header">
          <span class="slide-ref">Part 4 · Lecture Notes ${meta.slides}</span>
          <h2>4. Orthogonal Complements, Subspace Projections &amp; Best Approximation</h2>
          <p>How every vector $\\mathbf{y} \\in \\mathbb{R}^n$ splits uniquely into a shadow $\\hat{\\mathbf{y}} \\in W$ and a perpendicular residual $\\mathbf{z} \\in W^\\perp$—yielding the <strong>closest point</strong> in $W$ and explaining the geometry of the Four Fundamental Subspaces.</p>
        </div>

        <!-- SECTION A: CORE DEFINITIONS & THEOREMS -->
        <div class="card">
          <h3>Orthogonal Complement $W^\\perp$, Projection Formulae &amp; Best Approximation</h3>
          <div class="grid-2">
            <div class="callout definition">
              <strong>Orthogonal Complement $W^\\perp$:</strong><br/>
              For a subspace $W \\subseteq \\mathbb{R}^n$, its orthogonal complement is
              $$W^\\perp = \\{\\mathbf{z} \\in \\mathbb{R}^n : \\langle \\mathbf{z}, \\mathbf{w} \\rangle = 0 \\text{ for all } \\mathbf{w} \\in W\\}.$$
              <ul style="margin:6px 0 0 18px; padding:0; line-height:1.55; font-size:0.9rem;">
                <li>$W^\\perp$ is always a subspace, $(W^\\perp)^\\perp = W$, and $W \\cap W^\\perp = \\{\\mathbf{0}\\}$.</li>
                <li><strong>Dimension Law:</strong> $\\dim(W) + \\dim(W^\\perp) = n$.</li>
              </ul>
            </div>

            <div class="callout tip">
              <strong>Orthogonal Decomposition &amp; Best Approximation:</strong><br/>
              Every $\\mathbf{y} \\in \\mathbb{R}^n$ writes <em>uniquely</em> as $\\mathbf{y} = \\hat{\\mathbf{y}} + \\mathbf{z}$ with $\\hat{\\mathbf{y}} \\in W$ and $\\mathbf{z} \\in W^\\perp$.<br/>
              For every other $\\mathbf{w} \\in W$, Pythagoras gives:
              $$\\|\\mathbf{y} - \\mathbf{w}\\|^2 = \\|\\mathbf{y} - \\hat{\\mathbf{y}}\\|^2 + \\|\\hat{\\mathbf{y}} - \\mathbf{w}\\|^2 \\implies \\|\\mathbf{y} - \\hat{\\mathbf{y}}\\| < \\|\\mathbf{y} - \\mathbf{w}\\|.$$
            </div>
          </div>

          <div class="callout definition" style="margin-top:12px;">
            <strong>Three Ways to Compute $\\hat{\\mathbf{y}} = \\operatorname{proj}_W(\\mathbf{y})$ (Know When Each Applies!):</strong>
            <div style="display:grid; grid-template-columns:repeat(auto-fit, minmax(220px, 1fr)); gap:10px; margin-top:8px;">
              <div style="background:#FFFFFF; border:1px solid #E2E8F0; border-radius:8px; padding:9px 11px;">
                <strong style="color:#0F766E;">1. Orthonormal Basis $\\{\\mathbf{q}_i\\}$ ($U$)</strong>
                $$\\hat{\\mathbf{y}} = \\sum_{i=1}^k \\langle \\mathbf{y}, \\mathbf{q}_i \\rangle \\mathbf{q}_i = (U U^T)\\mathbf{y}$$
                <span style="font-size:0.82rem; color:var(--text-muted);">Projection matrix $P = U U^T$ ($P^T = P$, $P^2 = P$).</span>
              </div>
              <div style="background:#FFFFFF; border:1px solid #E2E8F0; border-radius:8px; padding:9px 11px;">
                <strong style="color:#1D4ED8;">2. Orthogonal Basis $\\{\\mathbf{u}_i\\}$</strong>
                $$\\hat{\\mathbf{y}} = \\sum_{i=1}^k \\frac{\\langle \\mathbf{y}, \\mathbf{u}_i \\rangle}{\\|\\mathbf{u}_i\\|^2}\\mathbf{u}_i$$
                <span style="font-size:0.82rem; color:var(--text-muted);"><em>Warning:</em> Fails if $\\{\\mathbf{u}_i\\}$ is not orthogonal!</span>
              </div>
              <div style="background:#FFFFFF; border:1px solid #E2E8F0; border-radius:8px; padding:9px 11px;">
                <strong style="color:#9333EA;">3. General Basis (Cols of $A$)</strong>
                $$\\hat{\\mathbf{y}} = A(A^T A)^{-1}A^T \\mathbf{y}$$
                <span style="font-size:0.82rem; color:var(--text-muted);">Works for ANY linearly independent basis (preview of Least Squares!).</span>
              </div>
            </div>
          </div>
        </div>

        <!-- SECTION B: INTERACTIVE LAB 1 — 2D LINE W vs W^perp + TANGENT CIRCLE -->
        <div class="card">
          <h3>Interactive Lab 1: 2D Decomposition $\\mathbf{y} = \\hat{\\mathbf{y}} + \\mathbf{z}$ &amp; Tangent-Circle Proof</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            Drag <strong style="color:#1D4ED8;">$\\mathbf{y}$ (blue)</strong> anywhere and drag the competitor <strong style="color:#DC2626;">$\\mathbf{w}$ (red dot)</strong> along the teal line $W$. The dashed circle centered at $\\mathbf{y}$ kisses $W$ tangentially at $\\hat{\\mathbf{y}}$—proving visually that every other $\\mathbf{w} \\in W$ lies strictly farther away!
          </p>

          <div class="m4-controls">
            <label style="font-size:0.86rem; font-weight:700;">Line $W$ Angle $\\alpha$: <span id="m4AlphaLbl" style="color:var(--polyu-red);">26.6°</span></label>
            <input type="range" id="m4AlphaSlider" min="-90" max="90" step="0.5" value="26.57" style="width:165px;" />
            <button class="m4-chip" data-alpha="0">Horizontal ($0^\\circ$)</button>
            <button class="m4-chip" data-alpha="26.57">Slope $\\frac{1}{2}$ ($\\mathbf{u}=(2,1)^T$)</button>
            <button class="m4-chip" data-alpha="45">$45^\\circ$ ($y=x$)</button>
            <button class="m4-chip" id="m4SnapW2D">Snap $\\mathbf{w} \\to \\hat{\\mathbf{y}}$</button>
          </div>

          <div class="grid-2" style="align-items:start;">
            <div class="m4-canvas-wrap">
              <canvas id="m4Canvas2D" class="m4-canvas"></canvas>
            </div>

            <div>
              <div class="m4-stat-grid">
                <div class="m4-stat">
                  <div class="m4-stat-label">Projection $\\hat{\\mathbf{y}} \\in W$</div>
                  <div class="m4-stat-val" id="m4Yhat2D" style="color:#0F766E;"></div>
                </div>
                <div class="m4-stat">
                  <div class="m4-stat-label">Residual $\\mathbf{z} = \\mathbf{y} - \\hat{\\mathbf{y}} \\in W^\\perp$</div>
                  <div class="m4-stat-val" id="m4Z2D" style="color:#B45309;"></div>
                </div>
                <div class="m4-stat">
                  <div class="m4-stat-label">Min Distance $\\|\\mathbf{y} - \\hat{\\mathbf{y}}\\|$</div>
                  <div class="m4-stat-val" id="m4DistMin2D" style="color:#15803D;"></div>
                </div>
                <div class="m4-stat">
                  <div class="m4-stat-label">Competitor $\\|\\mathbf{y} - \\mathbf{w}\\|$</div>
                  <div class="m4-stat-val" id="m4DistW2D" style="color:#DC2626;"></div>
                </div>
              </div>

              <div class="callout tip" style="margin-top:10px;" id="m4Pyth2D"></div>

              <div style="display:flex; gap:16px; flex-wrap:wrap; margin-top:10px; align-items:flex-start;">
                <div>
                  <div class="m4-stat-label" style="margin-bottom:4px;">$P_W = \\mathbf{q}_1\\mathbf{q}_1^T$</div>
                  <div id="m4PW2D"></div>
                </div>
                <div>
                  <div class="m4-stat-label" style="margin-bottom:4px;">$P_{W^\\perp} = I - P_W$</div>
                  <div id="m4PWp2D"></div>
                </div>
              </div>
              <div style="font-size:0.82rem; color:var(--text-muted); margin-top:6px;">
                Check: $P_W + P_{W^\\perp} = I_2$, $P_W^2 = P_W$, $P_W^T = P_W$, and $P_W P_{W^\\perp} = 0$.
              </div>
            </div>
          </div>
        </div>

        <!-- SECTION C: INTERACTIVE LAB 2 — 3D PLANE PROJECTION & BEST APPROXIMATION -->
        <div class="card">
          <h3>Interactive Lab 2: 3D Projection onto a Plane $W \\subset \\mathbb{R}^3$ &amp; Projection Matrix $P = U U^T$</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            Choose a plane $W$ in $\\mathbb{R}^3$, move the target vector $\\mathbf{y}$, and slide the competitor point $\\mathbf{w} = \\hat{\\mathbf{y}} + s\\mathbf{q}_1 + t\\mathbf{q}_2$ around inside the plane $W$ to watch the right triangle $(\\mathbf{y}, \\hat{\\mathbf{y}}, \\mathbf{w})$ in 3D!
          </p>

          <div class="m4-controls">
            <button class="m4-chip active" data-plane="zeroSum">Plane $x_1+x_2+x_3=0$</button>
            <button class="m4-chip" data-plane="exam">Exam Plane $\\operatorname{span}\\{(1,1,0)^T, (0,1,1)^T\\}$</button>
            <button class="m4-chip" data-plane="tilt">Tilted Plane $x_3 = 0.4x_1 - 0.3x_2$</button>
            <button class="m4-chip" data-plane="xy">$xy$-Plane ($x_3 = 0$)</button>
            <button class="m4-chip" id="m4SnapW3D">Snap $\\mathbf{w} \\to \\hat{\\mathbf{y}}$ ($s=t=0$)</button>
          </div>

          <div class="grid-2" style="align-items:start;">
            <div id="m4Plotly3D" class="m4-plotly"></div>

            <div>
              <div style="background:#F8FAFC; border:1px solid var(--border); border-radius:9px; padding:10px 12px; margin-bottom:10px;">
                <div class="m4-stat-label" style="margin-bottom:6px; color:#1D4ED8;">Target Vector $\\mathbf{y} = (y_1, y_2, y_3)^T$</div>
                <div class="m4-slider-row"><span>$y_1$</span><input type="range" id="m4Y1" min="-3.5" max="3.5" step="0.5" value="3" /><span id="m4Y1V">3</span></div>
                <div class="m4-slider-row"><span>$y_2$</span><input type="range" id="m4Y2" min="-3.5" max="3.5" step="0.5" value="1" /><span id="m4Y2V">1</span></div>
                <div class="m4-slider-row"><span>$y_3$</span><input type="range" id="m4Y3" min="-3.5" max="3.5" step="0.5" value="2" /><span id="m4Y3V">2</span></div>

                <div class="m4-stat-label" style="margin:10px 0 6px; color:#DC2626;">Competitor Offset $\\mathbf{w} = \\hat{\\mathbf{y}} + s\\mathbf{q}_1 + t\\mathbf{q}_2 \\in W$</div>
                <div class="m4-slider-row"><span>Shift $s$</span><input type="range" id="m4S3D" min="-2.5" max="2.5" step="0.25" value="1.25" /><span id="m4S3DV">1.25</span></div>
                <div class="m4-slider-row"><span>Shift $t$</span><input type="range" id="m4T3D" min="-2.5" max="2.5" step="0.25" value="-0.75" /><span id="m4T3DV">-0.75</span></div>
              </div>

              <div class="callout tip" id="m4Readout3D" style="margin-bottom:10px;"></div>

              <div style="display:flex; gap:14px; align-items:flex-start; flex-wrap:wrap;">
                <div>
                  <div class="m4-stat-label" style="margin-bottom:4px;">Projection Matrix $P = U U^T$</div>
                  <div id="m4PMat3D"></div>
                </div>
                <div style="flex:1; min-width:160px; font-size:0.84rem; color:#334155;" id="m4PCheck3D"></div>
              </div>
            </div>
          </div>

          <details class="code-details">
            <summary>Python (SymPy) Companion Code — Subspace Projection &amp; Distance to a Plane</summary>
            <pre><code>import sympy as sp

# Plane W = span{v1, v2} in R^3
v1 = sp.Matrix([1, 1, 0])
v2 = sp.Matrix([0, 1, 1])
y  = sp.Matrix([3, 1, 2])

# Method A: General basis matrix A = [v1, v2], P = A (A^T A)^(-1) A^T
A = sp.Matrix.hstack(v1, v2)
P = A * (A.T * A).inv() * A.T

y_hat = P * y          # Orthogonal projection onto W
z     = y - y_hat      # Residual in W^perp

print("P =\\n", P)
print("y_hat =", y_hat.T)
print("z     =", z.T)
print("Check z ⟂ v1, v2:", z.dot(v1), z.dot(v2))   # 0, 0
print("Check P^2 == P and P^T == P:", P*P == P, P.T == P)  # True, True
print("Shortest distance ||z|| =", z.norm())</code></pre>
          </details>
        </div>

        <!-- SECTION D: THE FOUR FUNDAMENTAL SUBSPACES EXPLORER -->
        <div class="card">
          <h3>The Four Fundamental Subspaces &amp; Orthogonal Splitting</h3>
          <div class="callout definition">
            <strong>Fundamental Theorem of Linear Algebra (Orthogonality Part):</strong><br/>
            For any $m \\times n$ matrix $A$, the domain $\\mathbb{R}^n$ and codomain $\\mathbb{R}^m$ each split into a pair of orthogonal complements:
            $$(\\operatorname{Row} A)^\\perp = \\operatorname{Nul}(A) \\text{ in } \\mathbb{R}^n, \\qquad (\\operatorname{Col} A)^\\perp = \\operatorname{Nul}(A^T) \\text{ in } \\mathbb{R}^m.$$
            Every input $\\mathbf{x} \\in \\mathbb{R}^n$ decomposes uniquely as $\\mathbf{x} = \\mathbf{x}_{\\text{row}} + \\mathbf{x}_{\\text{nul}}$, and $A$ kills the null component: $A\\mathbf{x} = A\\mathbf{x}_{\\text{row}} + \\mathbf{0} = A\\mathbf{x}_{\\text{row}}$.
          </div>

          <div class="m4-controls" style="margin-top:12px;">
            <strong style="font-size:0.86rem;">Select Matrix $A$:</strong>
            <button class="m4-chip active" data-fs="0">$A_1$ ($2\\times 3$, rank $1$)</button>
            <button class="m4-chip" data-fs="1">$A_2$ ($2\\times 3$, rank $2$)</button>
            <button class="m4-chip" data-fs="2">$A_3$ ($3\\times 3$, rank $2$)</button>
          </div>

          <div style="display:flex; gap:18px; flex-wrap:wrap; align-items:center; background:#FFFFFF; border:1px solid var(--border); border-radius:10px; padding:12px 14px;">
            <div id="m4FsMatTex"></div>
            <div style="flex:1; min-width:220px; font-size:0.88rem;" id="m4FsSummary"></div>
          </div>

          <div class="m4-fs-grid" id="m4FsBoxes"></div>

          <!-- Interactive x = x_row + x_nul tester -->
          <div style="margin-top:14px; background:#F8FAFC; border:1px solid var(--border); border-radius:10px; padding:12px 14px;">
            <div style="display:flex; flex-wrap:wrap; align-items:center; gap:10px; margin-bottom:10px;">
              <strong>Test the Domain Split $\\mathbf{x} = \\mathbf{x}_{\\text{row}} + \\mathbf{x}_{\\text{nul}}$:</strong>
              <span>$\\mathbf{x} = ($</span>
              <input type="number" class="m4-num-input" id="m4X1" value="3" step="1" />
              <span>,</span>
              <input type="number" class="m4-num-input" id="m4X2" value="1" step="1" />
              <span>,</span>
              <input type="number" class="m4-num-input" id="m4X3" value="4" step="1" />
              <span>$)^T$</span>
            </div>
            <div id="m4SplitResult" style="font-size:0.9rem; line-height:1.65;"></div>
          </div>
        </div>

        <!-- SECTION E: QUICK CONCEPT CHECK -->
        <div class="card">
          <h3>Quick Concept Check (Exam Pitfall Radar)</h3>

          <div class="m4-quiz-item" data-q="1" data-ans="C">
            <strong>Q1.</strong> Suppose $\\mathbf{v}_1 = \\begin{bmatrix}1\\\\0\\end{bmatrix}$ and $\\mathbf{v}_2 = \\begin{bmatrix}1\\\\1\\end{bmatrix}$ form a basis of $W = \\mathbb{R}^2$. To project $\\mathbf{y} = \\begin{bmatrix}1\\\\0\\end{bmatrix}$ onto $W$, a student writes
            $$\\hat{\\mathbf{y}} = \\frac{\\langle \\mathbf{y}, \\mathbf{v}_1\\rangle}{\\|\\mathbf{v}_1\\|^2}\\mathbf{v}_1 + \\frac{\\langle \\mathbf{y}, \\mathbf{v}_2\\rangle}{\\|\\mathbf{v}_2\\|^2}\\mathbf{v}_2 = \\begin{bmatrix}1\\\\0\\end{bmatrix} + \\frac{1}{2}\\begin{bmatrix}1\\\\1\\end{bmatrix} = \\begin{bmatrix}1.5\\\\0.5\\end{bmatrix} \\ne \\mathbf{y}.$$
            What went wrong?
            <div class="m4-opt-group">
              <button class="m4-opt-btn" data-opt="A">They forgot to normalize $\\mathbf{y}$ first.</button>
              <button class="m4-opt-btn" data-opt="B">You cannot project a vector onto $\\mathbb{R}^2$.</button>
              <button class="m4-opt-btn" data-opt="C">The basis $\\{\\mathbf{v}_1, \\mathbf{v}_2\\}$ is NOT orthogonal, so the sum formula is invalid.</button>
              <button class="m4-opt-btn" data-opt="D">They divided by $\\|\\mathbf{v}_2\\|^2$ instead of $\\|\\mathbf{v}_2\\|$.</button>
            </div>
            <div class="m4-fb" id="m4Fb1"></div>
          </div>

          <div class="m4-quiz-item" data-q="2" data-ans="B">
            <strong>Q2.</strong> Let $W$ be a $2$-dimensional subspace of $\\mathbb{R}^5$. What is $\\dim(W^\\perp)$?
            <div class="m4-opt-group">
              <button class="m4-opt-btn" data-opt="A">$2$</button>
              <button class="m4-opt-btn" data-opt="B">$3$</button>
              <button class="m4-opt-btn" data-opt="C">$5$</button>
              <button class="m4-opt-btn" data-opt="D">Depends on the choice of basis for $W$.</button>
            </div>
            <div class="m4-fb" id="m4Fb2"></div>
          </div>

          <div class="m4-quiz-item" data-q="3" data-ans="D">
            <strong>Q3.</strong> For any $m \\times n$ matrix $A$, which fundamental subspace is ALWAYS equal to $(\\operatorname{Col} A)^\\perp$?
            <div class="m4-opt-group">
              <button class="m4-opt-btn" data-opt="A">$\\operatorname{Nul}(A)$</button>
              <button class="m4-opt-btn" data-opt="B">$\\operatorname{Row}(A)$</button>
              <button class="m4-opt-btn" data-opt="C">$\\operatorname{Col}(A^T)$</button>
              <button class="m4-opt-btn" data-opt="D">$\\operatorname{Nul}(A^T)$</button>
            </div>
            <div class="m4-fb" id="m4Fb3"></div>
          </div>

          <div class="m4-quiz-item" data-q="4" data-ans="C">
            <strong>Q4.</strong> Suppose $\\mathbf{y} \\in \\mathbb{R}^4$ satisfies $\\|\\mathbf{y}\\| = 5$, and its orthogonal projection onto a subspace $W$ satisfies $\\|\\operatorname{proj}_W(\\mathbf{y})\\| = 3$. What is the <strong>shortest distance</strong> from $\\mathbf{y}$ to $W$?
            <div class="m4-opt-group">
              <button class="m4-opt-btn" data-opt="A">$2$</button>
              <button class="m4-opt-btn" data-opt="B">$\\sqrt{34}$</button>
              <button class="m4-opt-btn" data-opt="C">$4$</button>
              <button class="m4-opt-btn" data-opt="D">$16$</button>
            </div>
            <div class="m4-fb" id="m4Fb4"></div>
          </div>
        </div>
      `;
    },

    mount(container, { renderMath }) {
      // =====================================================================
      // LAB 1: 2D LINE W vs W^perp + TANGENT CIRCLE PROOF
      // =====================================================================
      const c2 = container.querySelector('#m4Canvas2D');
      const ctx2 = c2.getContext('2d');
      const alphaSlider = container.querySelector('#m4AlphaSlider');
      const alphaLbl = container.querySelector('#m4AlphaLbl');

      let alphaDeg = 26.57;
      let y2 = [1.5, 3.25];
      let wOffset2 = -1.6; // w = yhat + wOffset2 * q1
      let dragTarget2 = null; // 'y' or 'w'
      const range2 = 5.0;

      function layout2() {
        const r = c2.getBoundingClientRect();
        const w = r.width || 450, h = r.height || 370;
        const sc = Math.min(w, h) / (2 * range2);
        return { w, h, cx: w / 2, cy: h / 2, sc };
      }
      function toS2(p, L) { return [L.cx + p[0] * L.sc, L.cy - p[1] * L.sc]; }
      function toW2(sx, sy, L) {
        return [
          Math.max(-4.5, Math.min(4.5, (sx - L.cx) / L.sc)),
          Math.max(-4.5, Math.min(4.5, (L.cy - sy) / L.sc))
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
        const L = layout2();
        c2.width = Math.round(L.w * dpr);
        c2.height = Math.round(L.h * dpr);
        ctx2.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx2.clearRect(0, 0, L.w, L.h);

        const rad = (alphaDeg * Math.PI) / 180;
        const q1 = [Math.cos(rad), Math.sin(rad)];
        const q2 = [-Math.sin(rad), Math.cos(rad)];

        const yhat = scale(dot(y2, q1), q1);
        const z = sub(y2, yhat);
        const wPt = add(yhat, scale(wOffset2, q1));
        const dMin = norm(z);
        const dW = norm(sub(y2, wPt));

        // Faint coordinate grid
        ctx2.strokeStyle = '#F1F5F9'; ctx2.lineWidth = 1;
        for (let g = -5; g <= 5; g++) {
          const [gx] = toS2([g, 0], L), [, gy] = toS2([0, g], L);
          ctx2.beginPath(); ctx2.moveTo(gx, 0); ctx2.lineTo(gx, L.h); ctx2.stroke();
          ctx2.beginPath(); ctx2.moveTo(0, gy); ctx2.lineTo(L.w, gy); ctx2.stroke();
        }
        ctx2.strokeStyle = '#E2E8F0'; ctx2.lineWidth = 1.2;
        ctx2.beginPath(); ctx2.moveTo(0, L.cy); ctx2.lineTo(L.w, L.cy); ctx2.stroke();
        ctx2.beginPath(); ctx2.moveTo(L.cx, 0); ctx2.lineTo(L.cx, L.h); ctx2.stroke();

        // Subspace W (teal line) and W^perp (amber line)
        const wA = toS2(scale(-8, q1), L), wB = toS2(scale(8, q1), L);
        ctx2.save();
        ctx2.strokeStyle = 'rgba(15, 118, 110, 0.55)'; ctx2.lineWidth = 2.6;
        ctx2.beginPath(); ctx2.moveTo(wA[0], wA[1]); ctx2.lineTo(wB[0], wB[1]); ctx2.stroke();

        const pA = toS2(scale(-8, q2), L), pB = toS2(scale(8, q2), L);
        ctx2.strokeStyle = 'rgba(180, 83, 9, 0.42)'; ctx2.lineWidth = 1.8;
        ctx2.setLineDash([6, 4]);
        ctx2.beginPath(); ctx2.moveTo(pA[0], pA[1]); ctx2.lineTo(pB[0], pB[1]); ctx2.stroke();
        ctx2.restore();

        // Label W and W^perp
        const lblW = toS2(scale(3.9, q1), L);
        const lblWp = toS2(scale(3.6, q2), L);
        ctx2.font = 'bold 13px sans-serif';
        ctx2.fillStyle = '#0F766E'; ctx2.fillText('W', lblW[0] + 5, lblW[1] - 6);
        ctx2.fillStyle = '#B45309'; ctx2.fillText('W⊥', lblWp[0] + 5, lblWp[1] - 6);

        const [yx, yy] = toS2(y2, L);
        const [hx, hy] = toS2(yhat, L);
        const [zx, zy] = toS2(z, L);
        const [wx, wy] = toS2(wPt, L);

        // Tangent Circle centered at y with radius ||y - yhat||
        if (dMin > 0.05) {
          ctx2.save();
          ctx2.strokeStyle = 'rgba(21, 128, 61, 0.45)';
          ctx2.fillStyle = 'rgba(21, 128, 61, 0.06)';
          ctx2.lineWidth = 1.8;
          ctx2.setLineDash([5, 4]);
          ctx2.beginPath();
          ctx2.arc(yx, yy, dMin * L.sc, 0, Math.PI * 2);
          ctx2.fill();
          ctx2.stroke();
          ctx2.restore();
        }

        // Right-angle square at yhat
        if (dMin > 0.15) {
          const sq = 0.28;
          const dirZ = scale(1 / dMin, z);
          const dirW = wOffset2 >= 0 ? q1 : scale(-1, q1);
          const p1 = toS2(add(yhat, scale(sq, dirW)), L);
          const p2 = toS2(add(yhat, add(scale(sq, dirW), scale(sq, dirZ))), L);
          const p3 = toS2(add(yhat, scale(sq, dirZ)), L);
          ctx2.save();
          ctx2.strokeStyle = '#64748B'; ctx2.lineWidth = 1.5;
          ctx2.beginPath();
          ctx2.moveTo(p1[0], p1[1]); ctx2.lineTo(p2[0], p2[1]); ctx2.lineTo(p3[0], p3[1]);
          ctx2.stroke();
          ctx2.restore();
        }

        // Segment from w to y (competitor error)
        if (Math.abs(wOffset2) > 0.04) {
          ctx2.save();
          ctx2.strokeStyle = '#DC2626'; ctx2.lineWidth = 2.0; ctx2.setLineDash([5, 4]);
          ctx2.beginPath(); ctx2.moveTo(wx, wy); ctx2.lineTo(yx, yy); ctx2.stroke();
          ctx2.restore();
        }

        // Residual z at yhat -> y (amber dashed) and at origin -> z (amber solid)
        ctx2.save();
        ctx2.strokeStyle = '#D97706'; ctx2.lineWidth = 2.4; ctx2.setLineDash([4, 3]);
        ctx2.beginPath(); ctx2.moveTo(hx, hy); ctx2.lineTo(yx, yy); ctx2.stroke();
        ctx2.restore();
        drawArrow(ctx2, L.cx, L.cy, zx, zy, '#D97706', 2.4);

        // Vectors yhat (teal) and y (blue)
        drawArrow(ctx2, L.cx, L.cy, hx, hy, '#0F766E', 3.2);
        drawArrow(ctx2, L.cx, L.cy, yx, yy, '#1D4ED8', 3.2);

        // Draggable handle for w on line W
        ctx2.beginPath(); ctx2.arc(wx, wy, 7, 0, Math.PI * 2);
        ctx2.fillStyle = '#FEE2E2'; ctx2.fill();
        ctx2.lineWidth = 2.4; ctx2.strokeStyle = '#DC2626'; ctx2.stroke();
        ctx2.fillStyle = '#DC2626'; ctx2.font = 'bold 13px sans-serif';
        ctx2.fillText('w (drag)', wx + 9, wy + 15);

        // Draggable handle for y
        ctx2.beginPath(); ctx2.arc(yx, yy, 7, 0, Math.PI * 2);
        ctx2.fillStyle = '#FFFFFF'; ctx2.fill();
        ctx2.lineWidth = 2.6; ctx2.strokeStyle = '#1D4ED8'; ctx2.stroke();
        ctx2.fillStyle = '#1D4ED8'; ctx2.font = 'bold 14px sans-serif';
        ctx2.fillText('y', yx + 10, yy - 8);

        ctx2.fillStyle = '#0F766E'; ctx2.fillText('ŷ', hx + 8, hy - 8);
        ctx2.fillStyle = '#B45309'; ctx2.fillText('z', zx + 8, zy - 6);

        // Readouts
        alphaLbl.textContent = `${alphaDeg.toFixed(1)}°`;
        container.querySelector('#m4Yhat2D').textContent = fmtVec(yhat, 2);
        container.querySelector('#m4Z2D').textContent = fmtVec(z, 2);
        container.querySelector('#m4DistMin2D').textContent = dMin.toFixed(3);
        container.querySelector('#m4DistW2D').textContent = dW.toFixed(3);

        const diffHW = norm(sub(yhat, wPt));
        const pythEl = container.querySelector('#m4Pyth2D');
        const isMin = diffHW < 0.03;
        pythEl.innerHTML = isMin
          ? `<strong>Minimum Achieved!</strong> $\\mathbf{w} = \\hat{\\mathbf{y}}$, so $\\|\\mathbf{y} - \\mathbf{w}\\| = \\|\\mathbf{y} - \\hat{\\mathbf{y}}\\| = \\mathbf{${dMin.toFixed(3)}}$ (exact tangent point!).`
          : `<strong>Pythagoras Check:</strong> $\\|\\mathbf{y} - \\mathbf{w}\\|^2 = \\|\\mathbf{y} - \\hat{\\mathbf{y}}\\|^2 + \\|\\hat{\\mathbf{y}} - \\mathbf{w}\\|^2$<br/>
             $${(dW * dW).toFixed(2)} = ${(dMin * dMin).toFixed(2)} + ${(diffHW * diffHW).toFixed(2)} \\implies \\|\\mathbf{y}-\\mathbf{w}\\| = ${dW.toFixed(2)} > ${dMin.toFixed(2)}$`;
        renderMath(pythEl);

        const PW = outerProd(q1, q1);
        const PWp = outerProd(q2, q2);
        container.querySelector('#m4PW2D').innerHTML = renderMatHTML(PW, 2);
        container.querySelector('#m4PWp2D').innerHTML = renderMatHTML(PWp, 2);
      }

      function ptr2(e) {
        const r = c2.getBoundingClientRect();
        const cx = e.touches ? e.touches[0].clientX : e.clientX;
        const cy = e.touches ? e.touches[0].clientY : e.clientY;
        return [cx - r.left, cy - r.top];
      }
      function startDrag2(e) {
        const L = layout2(), [px, py] = ptr2(e);
        const rad = (alphaDeg * Math.PI) / 180;
        const q1 = [Math.cos(rad), Math.sin(rad)];
        const yhat = scale(dot(y2, q1), q1);
        const wPt = add(yhat, scale(wOffset2, q1));
        const [yx, yy] = toS2(y2, L), [wx, wy] = toS2(wPt, L);
        if (Math.hypot(px - yx, py - yy) < 24) {
          dragTarget2 = 'y'; e.preventDefault();
        } else if (Math.hypot(px - wx, py - wy) < 24) {
          dragTarget2 = 'w'; e.preventDefault();
        }
      }
      function moveDrag2(e) {
        if (!dragTarget2 || !c2.isConnected) return;
        if (e.touches) e.preventDefault();
        const L = layout2(), [px, py] = ptr2(e), pt = toW2(px, py, L);
        if (dragTarget2 === 'y') {
          y2 = [Math.round(pt[0] * 4) / 4, Math.round(pt[1] * 4) / 4];
        } else if (dragTarget2 === 'w') {
          const rad = (alphaDeg * Math.PI) / 180;
          const q1 = [Math.cos(rad), Math.sin(rad)];
          const yhat = scale(dot(y2, q1), q1);
          const off = dot(sub(pt, yhat), q1);
          wOffset2 = Math.max(-4, Math.min(4, Math.round(off * 10) / 10));
        }
        drawLab2D();
      }
      c2.addEventListener('mousedown', startDrag2);
      c2.addEventListener('touchstart', startDrag2, { passive: false });
      window.addEventListener('mousemove', moveDrag2);
      window.addEventListener('touchmove', moveDrag2, { passive: false });
      window.addEventListener('mouseup', () => { dragTarget2 = null; });
      window.addEventListener('touchend', () => { dragTarget2 = null; });

      alphaSlider.addEventListener('input', () => {
        alphaDeg = parseFloat(alphaSlider.value);
        drawLab2D();
      });
      container.querySelectorAll('[data-alpha]').forEach(btn => {
        btn.addEventListener('click', () => {
          alphaDeg = parseFloat(btn.dataset.alpha);
          alphaSlider.value = alphaDeg;
          drawLab2D();
        });
      });
      container.querySelector('#m4SnapW2D').addEventListener('click', () => {
        wOffset2 = 0;
        drawLab2D();
      });

      // =====================================================================
      // LAB 2: 3D PLANE PROJECTION & BEST APPROXIMATION (PLOTLY)
      // =====================================================================
      const planePresets = {
        zeroSum: { name: 'x₁ + x₂ + x₃ = 0', v1: [1, -1, 0], v2: [1, 0, -1] },
        exam:    { name: 'span{(1,1,0)ᵀ, (0,1,1)ᵀ}', v1: [1, 1, 0], v2: [0, 1, 1] },
        tilt:    { name: 'x₃ = 0.4x₁ - 0.3x₂', v1: [1, 0, 0.4], v2: [0, 1, -0.3] },
        xy:      { name: 'xy-plane (x₃ = 0)', v1: [1, 0, 0], v2: [0, 1, 0] }
      };
      let activePlaneKey = 'zeroSum';
      const p3dEl = container.querySelector('#m4Plotly3D');

      function updateLab3D() {
        const y = [
          parseFloat(container.querySelector('#m4Y1').value),
          parseFloat(container.querySelector('#m4Y2').value),
          parseFloat(container.querySelector('#m4Y3').value)
        ];
        const s = parseFloat(container.querySelector('#m4S3D').value);
        const t = parseFloat(container.querySelector('#m4T3D').value);
        container.querySelector('#m4Y1V').textContent = y[0];
        container.querySelector('#m4Y2V').textContent = y[1];
        container.querySelector('#m4Y3V').textContent = y[2];
        container.querySelector('#m4S3DV').textContent = s;
        container.querySelector('#m4T3DV').textContent = t;

        const preset = planePresets[activePlaneKey];
        const [q1, q2] = gramSchmidt([preset.v1, preset.v2]);
        const nVec = cross3(q1, q2); // unit normal spanning W^perp

        // P = q1 q1^T + q2 q2^T = U U^T
        const P = addMat(outerProd(q1, q1), outerProd(q2, q2));
        const yhat = matVec(P, y);
        const z = sub(y, yhat);
        const w = add(yhat, add(scale(s, q1), scale(t, q2)));

        const dMin = norm(z);
        const dW = norm(sub(y, w));
        const dBase = norm(sub(yhat, w));

        const readEl = container.querySelector('#m4Readout3D');
        readEl.innerHTML = `
          <strong>Orthogonal Decomposition in $\\mathbb{R}^3$:</strong><br/>
          $\\hat{\\mathbf{y}} = ${fmtVec(yhat, 2)} \\in W, \\quad \\mathbf{z} = ${fmtVec(z, 2)} \\in W^\\perp$ (check $\\langle \\hat{\\mathbf{y}}, \\mathbf{z} \\rangle = ${clean(dot(yhat, z), 3)}$)<br/>
          <span style="color:#15803D; font-weight:700;">Shortest $\\|\\mathbf{y}-\\hat{\\mathbf{y}}\\| = ${dMin.toFixed(3)}$</span> &nbsp;|&nbsp;
          <span style="color:#DC2626; font-weight:700;">Competitor $\\|\\mathbf{y}-\\mathbf{w}\\| = ${dW.toFixed(3)}$</span><br/>
          <span style="font-size:0.84rem;">Pythagoras: $${(dW * dW).toFixed(2)} = ${(dMin * dMin).toFixed(2)} + ${(dBase * dBase).toFixed(2)}$ ✓</span>
        `;
        renderMath(readEl);

        container.querySelector('#m4PMat3D').innerHTML = renderMatHTML(P, 2);
        const P2 = matMul(P, P);
        let maxErr = 0;
        for (let i = 0; i < 3; i++)
          for (let j = 0; j < 3; j++) maxErr = Math.max(maxErr, Math.abs(P2[i][j] - P[i][j]));

        const chkEl = container.querySelector('#m4PCheck3D');
        chkEl.innerHTML = `
          <strong>Matrix Properties:</strong><br/>
          • Symmetric ($P^T = P$): <strong style="color:#15803D;">✓</strong><br/>
          • Idempotent ($\\|P^2 - P\\|_{\\max} = ${maxErr.toFixed(4)}$): <strong style="color:#15803D;">✓</strong><br/>
          • Unit normal spanning $W^\\perp$:<br/>$\\mathbf{n} = ${fmtVec(nVec, 2)}$
        `;
        renderMath(chkEl);

        if (typeof Plotly === 'undefined') {
          p3dEl.innerHTML = '<div style="padding:20px; color:var(--text-muted);">3D plot requires Plotly.js.</div>';
          return;
        }

        // Plane surface grid
        const gridVals = [-3.6, 0, 3.6];
        const sx = [], sy = [], sz = [];
        for (const u of gridVals) {
          const rx = [], ry = [], rz = [];
          for (const v of gridVals) {
            const p = add(scale(u, q1), scale(v, q2));
            rx.push(p[0]); ry.push(p[1]); rz.push(p[2]);
          }
          sx.push(rx); sy.push(ry); sz.push(rz);
        }

        const traces = [
          // Plane W
          {
            type: 'surface', name: 'Plane W', showscale: false, opacity: 0.34,
            x: sx, y: sy, z: sz,
            colorscale: [[0, '#99F6E4'], [1, '#0F766E']],
            hoverinfo: 'skip'
          },
          // W^perp normal line through origin
          {
            type: 'scatter3d', mode: 'lines', name: 'W⊥ (normal line)',
            x: [-3.2 * nVec[0], 3.2 * nVec[0]],
            y: [-3.2 * nVec[1], 3.2 * nVec[1]],
            z: [-3.2 * nVec[2], 3.2 * nVec[2]],
            line: { color: '#9333EA', width: 3, dash: 'dot' }
          },
          // Vector y (blue)
          {
            type: 'scatter3d', mode: 'lines+markers+text', name: 'y (target)',
            x: [0, y[0]], y: [0, y[1]], z: [0, y[2]],
            text: ['', 'y'], textposition: 'top center', textfont: { color: '#1D4ED8', size: 13 },
            line: { color: '#1D4ED8', width: 6 }, marker: { size: [3, 6], color: '#1D4ED8' }
          },
          // Projection yhat (teal)
          {
            type: 'scatter3d', mode: 'lines+markers+text', name: 'ŷ = proj_W(y)',
            x: [0, yhat[0]], y: [0, yhat[1]], z: [0, yhat[2]],
            text: ['', 'ŷ'], textposition: 'bottom center', textfont: { color: '#0F766E', size: 13 },
            line: { color: '#0F766E', width: 6 }, marker: { size: [3, 6], color: '#0F766E' }
          },
          // Residual z = y - yhat (amber)
          {
            type: 'scatter3d', mode: 'lines', name: 'z = y - ŷ (⊥ W)',
            x: [yhat[0], y[0]], y: [yhat[1], y[1]], z: [yhat[2], y[2]],
            line: { color: '#D97706', width: 5, dash: 'dash' }
          },
          // Competitor w and error segment y - w (red)
          {
            type: 'scatter3d', mode: 'lines+markers+text', name: 'y - w (competitor)',
            x: [w[0], y[0]], y: [w[1], y[1]], z: [w[2], y[2]],
            text: ['w', ''], textposition: 'bottom center', textfont: { color: '#DC2626', size: 13 },
            line: { color: '#DC2626', width: 4, dash: 'dash' },
            marker: { size: [6, 1], color: '#DC2626' }
          },
          // Base leg yhat -> w on the plane
          {
            type: 'scatter3d', mode: 'lines', showlegend: false,
            x: [yhat[0], w[0]], y: [yhat[1], w[1]], z: [yhat[2], w[2]],
            line: { color: '#475569', width: 3, dash: 'dot' }
          }
        ];

        Plotly.react(p3dEl, traces, {
          margin: { l: 0, r: 0, t: 8, b: 0 },
          legend: { orientation: 'h', y: 0.99, x: 0.01, font: { size: 10.5 } },
          scene: {
            aspectmode: 'cube',
            xaxis: { range: [-4.2, 4.2], title: 'x₁', zerolinecolor: '#94A3B8' },
            yaxis: { range: [-4.2, 4.2], title: 'x₂', zerolinecolor: '#94A3B8' },
            zaxis: { range: [-4.2, 4.2], title: 'x₃', zerolinecolor: '#94A3B8' },
            camera: { eye: { x: 1.55, y: 1.45, z: 1.15 } }
          }
        }, { responsive: true, displayModeBar: false });
      }

      ['#m4Y1', '#m4Y2', '#m4Y3', '#m4S3D', '#m4T3D'].forEach(id => {
        container.querySelector(id).addEventListener('input', updateLab3D);
      });
      container.querySelectorAll('[data-plane]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-plane]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          activePlaneKey = btn.dataset.plane;
          updateLab3D();
        });
      });
      container.querySelector('#m4SnapW3D').addEventListener('click', () => {
        container.querySelector('#m4S3D').value = 0;
        container.querySelector('#m4T3D').value = 0;
        updateLab3D();
      });

      // =====================================================================
      // SECTION D: FOUR FUNDAMENTAL SUBSPACES EXPLORER
      // =====================================================================
      const fsPresets = [
        {
          name: 'A_1', m: 2, n: 3, r: 1,
          A: [[1, 2, 1], [2, 4, 2]],
          rowBasis: [[1, 2, 1]],
          nulBasis: [[-2, 1, 0], [-1, 0, 1]],
          colBasis: [[1, 2]],
          leftNulBasis: [[-2, 1]]
        },
        {
          name: 'A_2', m: 2, n: 3, r: 2,
          A: [[1, 0, 1], [0, 1, 1]],
          rowBasis: [[1, 0, 1], [0, 1, 1]],
          nulBasis: [[-1, -1, 1]],
          colBasis: [[1, 0], [0, 1]],
          leftNulBasis: []
        },
        {
          name: 'A_3', m: 3, n: 3, r: 2,
          A: [[1, 1, 0], [0, 1, 1], [1, 2, 1]],
          rowBasis: [[1, 1, 0], [0, 1, 1]],
          nulBasis: [[1, -1, 1]],
          colBasis: [[1, 0, 1], [1, 1, 2]],
          leftNulBasis: [[-1, -1, 1]]
        }
      ];
      let activeFsIdx = 0;

      function basisToLatex(vecs) {
        if (!vecs.length) return '\\{\\mathbf{0}\\}';
        return '\\operatorname{span}\\left\\{' +
          vecs.map(v => `\\begin{bmatrix}${v.join('\\\\')}\\end{bmatrix}`).join(',\\; ') +
          '\\right\\}';
      }

      function updateFsExplorer() {
        const P = fsPresets[activeFsIdx];
        const matLatex = `${P.name} = \\begin{bmatrix}` +
          P.A.map(row => row.join(' & ')).join('\\\\') + `\\end{bmatrix}`;

        const matEl = container.querySelector('#m4FsMatTex');
        matEl.innerHTML = `$$${matLatex}$$`;
        renderMath(matEl);

        const sumEl = container.querySelector('#m4FsSummary');
        sumEl.innerHTML = `
          <strong>Size:</strong> $${P.m} \\times ${P.n}$ &nbsp;|&nbsp; <strong>Rank:</strong> $r = ${P.r}$<br/>
          • In $\\mathbb{R}^{${P.n}}$: $\\dim(\\operatorname{Row} A) + \\dim(\\operatorname{Nul} A) = ${P.r} + ${P.n - P.r} = ${P.n}$<br/>
          • In $\\mathbb{R}^{${P.m}}$: $\\dim(\\operatorname{Col} A) + \\dim(\\operatorname{Nul} A^T) = ${P.r} + ${P.m - P.r} = ${P.m}$
        `;
        renderMath(sumEl);

        const boxesEl = container.querySelector('#m4FsBoxes');
        boxesEl.innerHTML = `
          <div class="m4-fs-box" style="border-top:3px solid #0F766E;">
            <h4><span style="color:#0F766E;">$\\operatorname{Row}(A) \\subseteq \\mathbb{R}^{${P.n}}$</span>
              <span class="m4-dim-pill" style="background:#CCFBF1; color:#115E59;">dim = ${P.r}</span></h4>
            <div>$${basisToLatex(P.rowBasis)}$</div>
          </div>
          <div class="m4-fs-box" style="border-top:3px solid #B45309;">
            <h4><span style="color:#B45309;">$\\operatorname{Nul}(A) = (\\operatorname{Row} A)^\\perp \\subseteq \\mathbb{R}^{${P.n}}$</span>
              <span class="m4-dim-pill" style="background:#FEF3C7; color:#92400E;">dim = ${P.n - P.r}</span></h4>
            <div>$${basisToLatex(P.nulBasis)}$</div>
          </div>
          <div class="m4-fs-box" style="border-top:3px solid #1D4ED8;">
            <h4><span style="color:#1D4ED8;">$\\operatorname{Col}(A) \\subseteq \\mathbb{R}^{${P.m}}$</span>
              <span class="m4-dim-pill" style="background:#DBEAFE; color:#1E40AF;">dim = ${P.r}</span></h4>
            <div>$${basisToLatex(P.colBasis)}$</div>
          </div>
          <div class="m4-fs-box" style="border-top:3px solid #9333EA;">
            <h4><span style="color:#9333EA;">$\\operatorname{Nul}(A^T) = (\\operatorname{Col} A)^\\perp \\subseteq \\mathbb{R}^{${P.m}}$</span>
              <span class="m4-dim-pill" style="background:#F3E8FF; color:#6B21A8;">dim = ${P.m - P.r}</span></h4>
            <div>$${basisToLatex(P.leftNulBasis)}$</div>
          </div>
        `;
        renderMath(boxesEl);
        updateSplitTest();
      }

      function updateSplitTest() {
        const P = fsPresets[activeFsIdx];
        const x = [
          parseFloat(container.querySelector('#m4X1').value) || 0,
          parseFloat(container.querySelector('#m4X2').value) || 0,
          parseFloat(container.querySelector('#m4X3').value) || 0
        ];
        // Project x onto Row(A) using Gram-Schmidt orthonormal basis of Row(A)
        const qRow = gramSchmidt(P.rowBasis);
        let xRow = [0, 0, 0];
        for (const q of qRow) xRow = add(xRow, scale(dot(x, q), q));
        const xNul = sub(x, xRow);

        const Ax = matVec(P.A, x);
        const AxRow = matVec(P.A, xRow);
        const AxNul = matVec(P.A, xNul);

        const resEl = container.querySelector('#m4SplitResult');
        resEl.innerHTML = `
          • $\\mathbf{x}_{\\text{row}} = ${fmtVec(xRow, 2)} \\in \\operatorname{Row}(A)$, &nbsp;
            $\\mathbf{x}_{\\text{nul}} = ${fmtVec(xNul, 2)} \\in \\operatorname{Nul}(A)$<br/>
          • Orthogonality check: $\\langle \\mathbf{x}_{\\text{row}}, \\mathbf{x}_{\\text{nul}} \\rangle = \\mathbf{${clean(dot(xRow, xNul), 3)}}$ ✓<br/>
          • Action of $A$: $A\\mathbf{x}_{\\text{nul}} = ${fmtVec(AxNul, 2)}$ (killed!), so
            <strong style="color:#1D4ED8;">$A\\mathbf{x} = ${fmtVec(Ax, 2)} = A\\mathbf{x}_{\\text{row}}$</strong> ✓
        `;
        renderMath(resEl);
      }

      container.querySelectorAll('[data-fs]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-fs]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          activeFsIdx = parseInt(btn.dataset.fs, 10);
          updateFsExplorer();
        });
      });
      ['#m4X1', '#m4X2', '#m4X3'].forEach(id => {
        container.querySelector(id).addEventListener('input', updateSplitTest);
      });

      // =====================================================================
      // SECTION E: QUICK CONCEPT CHECK
      // =====================================================================
      const explanations = {
        1: 'Because ⟨v₁, v₂⟩ = 1 ≠ 0, {v₁, v₂} is NOT orthogonal! The sum of 1D projections only produces the true orthogonal projection onto W when the basis vectors are mutually orthogonal.',
        2: 'By the Dimension Law for orthogonal complements, dim(W) + dim(W⊥) = n = 5, so dim(W⊥) = 5 - 2 = 3.',
        3: 'A vector z ∈ ℝᵐ is orthogonal to every column of A iff aᵢᵀz = 0 for all columns i, which is exactly Aᵀz = 0, i.e., z ∈ Nul(Aᵀ). (Nul(A) lives in ℝⁿ and is orthogonal to Row(A)!)',
        4: 'By the Orthogonal Decomposition Theorem, y = ŷ + z with ŷ ⊥ z. Pythagoras gives ‖y‖² = ‖ŷ‖² + ‖z‖², so ‖z‖ = √(5² - 3²) = √16 = 4 (NOT 5 - 3 = 2!).'
      };

      container.querySelectorAll('.m4-quiz-item').forEach(item => {
        const qNum = item.dataset.q;
        const ans = item.dataset.ans;
        const fb = container.querySelector(`#m4Fb${qNum}`);
        item.querySelectorAll('.m4-opt-btn').forEach(btn => {
          btn.addEventListener('click', () => {
            item.querySelectorAll('.m4-opt-btn').forEach(b => b.classList.remove('correct', 'wrong'));
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
      updateFsExplorer();
    }
  });
})();
