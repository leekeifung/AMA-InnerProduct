/**
 * AMA1751 Linear Algebra — Topic 6: Inner Product
 * Module 2: Orthogonality, Projections & Cauchy-Schwarz (Slides P. 5–10)
 */
(function () {
  // Local vector helpers (self-contained so we only depend on LA.dot/norm/add/scale)
  function sub(a, b) { return a.map((v, i) => v - b[i]); }
  function gcd(a, b) {
    a = Math.abs(Math.round(a));
    b = Math.abs(Math.round(b));
    while (b) { const t = b; b = a % b; a = t; }
    return a || 1;
  }
  function isNearInt(x) { return Math.abs(x - Math.round(x)) < 1e-9; }

  // Format num/den as exact reduced LaTeX fraction when both are integers
  function fmtFrac(num, den) {
    if (Math.abs(den) < 1e-12) return '\\text{undef}';
    if (isNearInt(num) && isNearInt(den)) {
      let n = Math.round(num), d = Math.round(den);
      if (d < 0) { n = -n; d = -d; }
      if (n === 0) return '0';
      const g = gcd(n, d);
      n /= g; d /= g;
      if (d === 1) return `${n}`;
      return n < 0 ? `-\\frac{${-n}}{${d}}` : `\\frac{${n}}{${d}}`;
    }
    return (num / den).toFixed(3);
  }

  function fmtVecExact(arr, den) {
    if (isNearInt(den) && arr.every(isNearInt)) {
      return arr.map(v => fmtFrac(v, den)).join(' \\\\ ');
    }
    return arr.map(v => (v / den).toFixed(3)).join(' \\\\ ');
  }

  AMA1751.registerModule(2, {
    render(meta) {
      return `
        <style>
          .m2-canvas-wrap {
            position: relative;
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 10px;
            overflow: hidden;
          }
          .m2-canvas {
            display: block;
            width: 100%;
            height: 370px;
            cursor: grab;
            touch-action: none;
          }
          .m2-canvas:active { cursor: grabbing; }
          .m2-controls {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            margin-bottom: 12px;
            align-items: center;
          }
          .m2-chip {
            padding: 6px 11px;
            border-radius: 999px;
            border: 1px solid var(--border);
            background: #F8FAFC;
            font-size: 0.82rem;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.15s;
          }
          .m2-chip:hover { background: var(--polyu-red-soft); border-color: #E8B4C4; color: var(--polyu-red-dark); }
          .m2-stat-grid {
            display: grid;
            grid-template-columns: repeat(2, minmax(0, 1fr));
            gap: 10px;
          }
          .m2-stat {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 10px 12px;
          }
          .m2-stat-label { font-size: 0.75rem; color: var(--text-muted); font-weight: 600; text-transform: uppercase; letter-spacing: 0.04em; }
          .m2-stat-val { font-size: 0.98rem; font-weight: 700; margin-top: 3px; color: #0F172A; }
          .m2-vec-input {
            width: 64px;
            padding: 6px 8px;
            border: 1px solid #CBD5E1;
            border-radius: 6px;
            text-align: center;
            font-weight: 600;
            font-size: 0.92rem;
          }
          .m2-vec-input:focus { outline: 2px solid var(--polyu-red); border-color: var(--polyu-red); }
          .m2-plotly {
            width: 100%;
            height: 390px;
            border: 1px solid var(--border);
            border-radius: 10px;
            background: #FFFFFF;
            overflow: hidden;
          }
          .m2-math-scroll { overflow-x: auto; padding: 12px 16px; background: #F8FAFC; border: 1px solid var(--border); border-radius: 10px; }
          .m2-meter-track {
            background: #E2E8F0;
            height: 12px;
            border-radius: 999px;
            overflow: hidden;
            margin: 8px 0 4px;
          }
          .m2-meter-fill {
            height: 100%;
            background: linear-gradient(90deg, #0F766E, #1D4ED8);
            transition: width 0.15s;
          }
          .m2-quiz-item {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 14px 16px;
            margin-bottom: 12px;
          }
          .m2-opt-group { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 10px; }
          .m2-opt-btn {
            padding: 7px 13px;
            border: 1px solid #CBD5E1;
            border-radius: 8px;
            background: #FFFFFF;
            cursor: pointer;
            font-weight: 600;
            font-size: 0.88rem;
          }
          .m2-opt-btn:hover { background: #F1F5F9; }
          .m2-opt-btn.correct { background: #DCFCE7; border-color: #16A34A; color: #166534; }
          .m2-opt-btn.wrong { background: #FEE2E2; border-color: #DC2626; color: #991B1B; }
          .m2-fb { margin-top: 8px; font-size: 0.88rem; font-weight: 500; }
        </style>

        <!-- HEADER -->
        <div class="module-header">
          <span class="slide-ref">Part 2 · Lecture Notes ${meta.slides}</span>
          <h2>2. Orthogonality, Projections &amp; Cauchy–Schwarz</h2>
          <p>Perpendicularity in $\\mathbb{R}^n$, decomposing any vector into a <strong>parallel projection</strong> and an <strong>orthogonal residual</strong>, and defining angles in any dimension via the Cauchy–Schwarz inequality.</p>
        </div>

        <!-- SECTION A: DEFINITIONS & PYTHAGOREAN THEOREM -->
        <div class="card">
          <h3>Orthogonal Vectors &amp; The Pythagorean Theorem in $\\mathbb{R}^n$</h3>
          <div class="grid-2">
            <div class="callout definition">
              <strong>Definition — Orthogonality:</strong><br/>
              Two vectors $\\mathbf{x}, \\mathbf{y} \\in \\mathbb{R}^n$ are <strong>orthogonal</strong> (written $\\mathbf{x} \\perp \\mathbf{y}$) if
              $$\\langle \\mathbf{x}, \\mathbf{y} \\rangle = \\mathbf{y}^T \\mathbf{x} = 0.$$
              <span style="font-size:0.88rem; color:var(--text-muted);">Note: The zero vector $\\mathbf{0}$ is orthogonal to <em>every</em> vector in $\\mathbb{R}^n$.</span>
            </div>

            <div class="callout tip">
              <strong>Theorem — Pythagorean Theorem in $\\mathbb{R}^n$:</strong><br/>
              Two vectors $\\mathbf{x}$ and $\\mathbf{y}$ are orthogonal <em>if and only if</em>
              $$\\|\\mathbf{x} + \\mathbf{y}\\|^2 = \\|\\mathbf{x}\\|^2 + \\|\\mathbf{y}\\|^2.$$
              <span style="font-size:0.86rem; color:var(--text-muted);"><strong>Why?</strong> Expanding $\\|\\mathbf{x}+\\mathbf{y}\\|^2 = \\langle \\mathbf{x}+\\mathbf{y}, \\mathbf{x}+\\mathbf{y}\\rangle = \\|\\mathbf{x}\\|^2 + 2\\langle \\mathbf{x},\\mathbf{y}\\rangle + \\|\\mathbf{y}\\|^2$, the cross term vanishes precisely when $\\langle \\mathbf{x},\\mathbf{y}\\rangle = 0$.</span>
            </div>
          </div>

          <div class="callout definition" style="margin-top:12px;">
            <strong>Orthogonal Projection of $\\mathbf{x}$ onto a nonzero vector $\\mathbf{y}$:</strong><br/>
            We can uniquely split $\\mathbf{x} = \\mathbf{p} + \\mathbf{w}$, where $\\mathbf{p} \\in \\operatorname{span}\\{\\mathbf{y}\\}$ (parallel to $\\mathbf{y}$) and $\\mathbf{w} \\perp \\mathbf{y}$:
            $$\\mathbf{p} = \\operatorname{proj}_{\\mathbf{y}}\\mathbf{x} = \\frac{\\langle \\mathbf{x}, \\mathbf{y} \\rangle}{\\langle \\mathbf{y}, \\mathbf{y} \\rangle}\\,\\mathbf{y} = \\frac{\\langle \\mathbf{x}, \\mathbf{y} \\rangle}{\\|\\mathbf{y}\\|^2}\\,\\mathbf{y}, \\qquad \\mathbf{w} = \\mathbf{x} - \\mathbf{p}.$$
            <span style="font-size:0.88rem; color:var(--text-muted);"><strong>Geometric meaning:</strong> $\\mathbf{p}$ is the <strong>closest point</strong> on the line $\\operatorname{span}\\{\\mathbf{y}\\}$ to $\\mathbf{x}$, and $\\|\\mathbf{w}\\| = \\|\\mathbf{x} - \\mathbf{p}\\|$ is the shortest distance from $\\mathbf{x}$ to that line.</span>
          </div>
        </div>

        <!-- SECTION B: INTERACTIVE LAB 1 — 2D PROJECTION & CLOSEST-POINT EXPLORER -->
        <div class="card">
          <h3>Interactive Lab 1: 2D Projection &amp; "Closest Point" Explorer</h3>
          <p style="margin-top:0; color:var(--text-muted);">Drag <strong style="color:#1D4ED8;">$\\mathbf{x}$ (blue)</strong> or <strong style="color:#DC2626;">$\\mathbf{y}$ (red)</strong>. Watch how <strong style="color:#0F766E;">$\\mathbf{p} = \\operatorname{proj}_{\\mathbf{y}}\\mathbf{x}$</strong> always sits at the foot of the perpendicular from $\\mathbf{x}$ onto $\\operatorname{span}\\{\\mathbf{y}\\}$.</p>

          <div class="m2-controls">
            <button class="m2-chip" data-preset="acute">Acute Projection</button>
            <button class="m2-chip" data-preset="obtuse">Obtuse ($\\mathbf{p}$ flips opposite $\\mathbf{y}$)</button>
            <button class="m2-chip" data-preset="ortho">Orthogonal ($\\mathbf{p}=\\mathbf{0}$)</button>
            <button class="m2-chip" data-preset="parallel">Parallel ($\\mathbf{w}=\\mathbf{0}$)</button>
            <label class="m2-chip" style="display:inline-flex; align-items:center; gap:5px;">
              <input type="checkbox" id="m2Snap" checked /> Snap to 0.5
            </label>
          </div>

          <div class="grid-2" style="align-items:start;">
            <div class="m2-canvas-wrap">
              <canvas id="m2Canvas" class="m2-canvas"></canvas>
            </div>

            <div>
              <div class="m2-stat-grid">
                <div class="m2-stat">
                  <div class="m2-stat-label">Projection $\\mathbf{p} = \\operatorname{proj}_{\\mathbf{y}}\\mathbf{x}$</div>
                  <div class="m2-stat-val" id="m2ProjVal" style="color:#0F766E;"></div>
                </div>
                <div class="m2-stat">
                  <div class="m2-stat-label">Residual $\\mathbf{w} = \\mathbf{x} - \\mathbf{p}$</div>
                  <div class="m2-stat-val" id="m2ResVal" style="color:#B45309;"></div>
                </div>
                <div class="m2-stat">
                  <div class="m2-stat-label">Orthogonality Check $\\langle \\mathbf{w}, \\mathbf{y} \\rangle$</div>
                  <div class="m2-stat-val" id="m2OrthoChk"></div>
                </div>
                <div class="m2-stat">
                  <div class="m2-stat-label">Pythagorean Check</div>
                  <div class="m2-stat-val" id="m2PythChk" style="font-size:0.86rem;"></div>
                </div>
              </div>

              <!-- Closest-Point Discovery Slider -->
              <div class="callout tip" style="margin-top:12px;">
                <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:6px;">
                  <strong>Why Projection = Best Approximation!</strong>
                  <label style="font-size:0.82rem; font-weight:600; cursor:pointer;">
                    <input type="checkbox" id="m2ShowCandidate" checked /> Show test point $t\\mathbf{y}$
                  </label>
                </div>
                <p style="margin:6px 0 4px; font-size:0.85rem; color:#334155;">
                  Slide $t$ to move the purple test point $t\\mathbf{y}$ along $\\operatorname{span}\\{\\mathbf{y}\\}$ and compare its error $\\|\\mathbf{x} - t\\mathbf{y}\\|$ against the projection error $\\|\\mathbf{x} - \\mathbf{p}\\|$:
                </p>
                <div style="display:flex; align-items:center; gap:10px;">
                  <input type="range" id="m2TSlider" min="-2.5" max="2.5" step="0.05" value="1.5" style="flex:1;" />
                  <button class="m2-chip" id="m2SnapOptimal" style="white-space:nowrap;">Jump to $t^*$</button>
                </div>
                <div id="m2ClosestInfo" style="margin-top:6px; font-size:0.86rem; font-weight:600;"></div>
              </div>
            </div>
          </div>
        </div>

        <!-- SECTION C: INTERACTIVE LAB 2 — 3D PLOTLY PROJECTION & EXAMPLE 6.1 -->
        <div class="card">
          <h3>Interactive Lab 2: 3D Vector Projection (Plotly) &amp; Worked Example</h3>
          <p style="margin-top:0; color:var(--text-muted);">Decompose $\\mathbf{x} \\in \\mathbb{R}^3$ into $\\mathbf{p} \\parallel \\mathbf{y}$ and $\\mathbf{w} \\perp \\mathbf{y}$. Rotate the 3D scene with your mouse/finger to verify the $90^\\circ$ angle from every viewpoint!</p>

          <div class="m2-controls">
            <button class="btn" id="m2PresetEx1">Example 6.1: $\\mathbf{x}=(1,2,3)^T, \\mathbf{y}=(1,1,1)^T$</button>
            <button class="btn" id="m2PresetEx2">Exact Integer $\\mathbf{p}$: $\\mathbf{x}=(4,1,2)^T, \\mathbf{y}=(1,2,2)^T$</button>
            <button class="btn" id="m2PresetEx3">Obtuse 3D: $\\mathbf{x}=(-2,3,1)^T, \\mathbf{y}=(2,1,-1)^T$</button>
          </div>

          <div class="m2-controls" style="gap:16px;">
            <div style="display:flex; align-items:center; gap:6px;">
              <strong style="color:#1D4ED8;">$\\mathbf{x}^T =$</strong>
              <input type="number" step="any" class="m2-vec-input m2-3d-x" value="1" />
              <input type="number" step="any" class="m2-vec-input m2-3d-x" value="2" />
              <input type="number" step="any" class="m2-vec-input m2-3d-x" value="3" />
            </div>
            <div style="display:flex; align-items:center; gap:6px;">
              <strong style="color:#DC2626;">$\\mathbf{y}^T =$</strong>
              <input type="number" step="any" class="m2-vec-input m2-3d-y" value="1" />
              <input type="number" step="any" class="m2-vec-input m2-3d-y" value="1" />
              <input type="number" step="any" class="m2-vec-input m2-3d-y" value="1" />
            </div>
          </div>

          <div class="grid-2" style="align-items:start;">
            <div id="m2Plotly3D" class="m2-plotly"></div>
            <div id="m2Calc3D" class="m2-math-scroll"></div>
          </div>

          <details class="code-details">
            <summary>Python (SymPy) Companion Code — Orthogonal Projection &amp; Verification</summary>
            <pre><code>import sympy as sp

x = sp.Matrix([1, 2, 3])
y = sp.Matrix([1, 1, 1])

# Projection of x onto y: p = (<x, y> / <y, y>) * y
p = (x.dot(y) / y.dot(y)) * y
w = x - p

print("Projection p =", p.T)          # Matrix([[2, 2, 2]])
print("Residual   w =", w.T)          # Matrix([[-1, 0, 1]])
print("<w, y>      =", w.dot(y))      # 0  (confirms w ⊥ y)
print("||p||^2 + ||w||^2 =", p.norm()**2 + w.norm()**2, "vs ||x||^2 =", x.norm()**2)</code></pre>
          </details>
        </div>

        <!-- SECTION D: CAUCHY-SCHWARZ INEQUALITY & ANGLE IN R^n -->
        <div class="card">
          <h3>Cauchy–Schwarz Inequality &amp; Angle Between Vectors in $\\mathbb{R}^n$</h3>
          <div class="grid-2">
            <div class="callout definition">
              <strong>Theorem — Cauchy–Schwarz Inequality:</strong><br/>
              For all $\\mathbf{x}, \\mathbf{y} \\in \\mathbb{R}^n$,
              $$|\\langle \\mathbf{x}, \\mathbf{y} \\rangle| \\le \\|\\mathbf{x}\\|\\,\\|\\mathbf{y}\\|.$$
              Equality holds <em>if and only if</em> $\\mathbf{x}$ and $\\mathbf{y}$ are linearly dependent (one is a scalar multiple of the other, or one is $\\mathbf{0}$).
            </div>

            <div class="callout tip">
              <strong>Consequence — Angle $\\theta$ in ANY Dimension $\\mathbb{R}^n$:</strong><br/>
              Dividing by $\\|\\mathbf{x}\\|\\|\\mathbf{y}\\|$ (for nonzero vectors) gives $-1 \\le \\frac{\\langle \\mathbf{x},\\mathbf{y}\\rangle}{\\|\\mathbf{x}\\|\\|\\mathbf{y}\\|} \\le 1$. Thus there is a unique angle $\\theta \\in [0, \\pi]$ satisfying:
              $$\\cos \\theta = \\frac{\\langle \\mathbf{x}, \\mathbf{y} \\rangle}{\\|\\mathbf{x}\\|\\,\\|\\mathbf{y}\\|}.$$
              <span style="font-size:0.84rem; color:var(--text-muted);">In AI &amp; data science, $\\cos\\theta$ is called <strong>Cosine Similarity</strong> (used to compare high-dimensional document/user embeddings!).</span>
            </div>
          </div>

          <p style="font-size:0.9rem; color:#334155; margin:10px 0;">
            <strong>2-Line Proof from Projection:</strong> Because $\\mathbf{w} = \\mathbf{x} - \\mathbf{p}$ has non-negative squared length,
            $0 \\le \\|\\mathbf{x} - \\mathbf{p}\\|^2 = \\|\\mathbf{x}\\|^2 - \\|\\mathbf{p}\\|^2 = \\|\\mathbf{x}\\|^2 - \\frac{\\langle \\mathbf{x},\\mathbf{y}\\rangle^2}{\\|\\mathbf{y}\\|^2} \\implies \\langle \\mathbf{x},\\mathbf{y}\\rangle^2 \\le \\|\\mathbf{x}\\|^2\\|\\mathbf{y}\\|^2$.
          </p>

          <!-- Interactive Cauchy-Schwarz & Cosine Similarity Tester in R^4 -->
          <div style="background:#F8FAFC; border:1px solid var(--border); border-radius:10px; padding:14px 16px; margin-top:10px;">
            <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px; margin-bottom:10px;">
              <strong>Interactive Cauchy–Schwarz &amp; Cosine Similarity Tester ($\\mathbb{R}^4$)</strong>
              <div style="display:flex; flex-wrap:wrap; gap:6px;">
                <button class="m2-chip" data-cs="parallel">Parallel (Equality!)</button>
                <button class="m2-chip" data-cs="anti">Anti-Parallel (Equality!)</button>
                <button class="m2-chip" data-cs="ortho">Orthogonal ($90^\\circ$)</button>
                <button class="m2-chip" data-cs="random">Random $\\mathbb{R}^4$</button>
              </div>
            </div>

            <div class="m2-controls" style="gap:14px; margin-bottom:10px;">
              <div style="display:flex; align-items:center; gap:5px; flex-wrap:wrap;">
                <strong style="color:#1D4ED8;">$\\mathbf{u}^T =$</strong>
                <input type="number" step="any" class="m2-vec-input m2-cs-u" value="1" />
                <input type="number" step="any" class="m2-vec-input m2-cs-u" value="2" />
                <input type="number" step="any" class="m2-vec-input m2-cs-u" value="0" />
                <input type="number" step="any" class="m2-vec-input m2-cs-u" value="-2" />
              </div>
              <div style="display:flex; align-items:center; gap:5px; flex-wrap:wrap;">
                <strong style="color:#DC2626;">$\\mathbf{v}^T =$</strong>
                <input type="number" step="any" class="m2-vec-input m2-cs-v" value="2" />
                <input type="number" step="any" class="m2-vec-input m2-cs-v" value="1" />
                <input type="number" step="any" class="m2-vec-input m2-cs-v" value="3" />
                <input type="number" step="any" class="m2-vec-input m2-cs-v" value="-1" />
              </div>
            </div>

            <div id="m2CsReadout" style="font-size:0.93rem;"></div>
            <div class="m2-meter-track"><div id="m2CsBar" class="m2-meter-fill" style="width:50%;"></div></div>
            <div style="font-size:0.79rem; color:var(--text-muted);">
              Bar shows ratio $|\\langle \\mathbf{u},\\mathbf{v}\\rangle| / (\\|\\mathbf{u}\\|\\|\\mathbf{v}\\|) = |\\cos\\theta| \\le 100\\%$. It reaches $100\\%$ only when $\\mathbf{u}$ and $\\mathbf{v}$ are parallel!
            </div>
          </div>
        </div>

        <!-- SECTION E: QUICK CONCEPT CHECK -->
        <div class="card">
          <h3>Quick Concept Check (Exam Pitfall Radar)</h3>

          <div class="m2-quiz-item" data-q="1" data-ans="C">
            <strong>Q1.</strong> For what value of $k$ are $\\mathbf{x} = \\begin{bmatrix} 2 \\\\ k \\\\ -3 \\end{bmatrix}$ and $\\mathbf{y} = \\begin{bmatrix} 4 \\\\ -2 \\\\ 2 \\end{bmatrix}$ orthogonal?
            <div class="m2-opt-group">
              <button class="m2-opt-btn" data-opt="A">$k = -1$</button>
              <button class="m2-opt-btn" data-opt="B">$k = 2$</button>
              <button class="m2-opt-btn" data-opt="C">$k = 1$</button>
              <button class="m2-opt-btn" data-opt="D">$k = 7$</button>
            </div>
            <div class="m2-fb" id="m2Fb1"></div>
          </div>

          <div class="m2-quiz-item" data-q="2" data-ans="B">
            <strong>Q2.</strong> Find the orthogonal projection $\\operatorname{proj}_{\\mathbf{y}}\\mathbf{x}$ of $\\mathbf{x} = \\begin{bmatrix} 5 \\\\ 0 \\end{bmatrix}$ onto $\\mathbf{y} = \\begin{bmatrix} 3 \\\\ 4 \\end{bmatrix}$.
            <div class="m2-opt-group">
              <button class="m2-opt-btn" data-opt="A">$3\\begin{bmatrix} 3 \\\\ 4 \\end{bmatrix} = \\begin{bmatrix} 9 \\\\ 12 \\end{bmatrix}$</button>
              <button class="m2-opt-btn" data-opt="B">$\\frac{3}{5}\\begin{bmatrix} 3 \\\\ 4 \\end{bmatrix} = \\begin{bmatrix} 9/5 \\\\ 12/5 \\end{bmatrix}$</button>
              <button class="m2-opt-btn" data-opt="C">$\\frac{3}{5}\\begin{bmatrix} 5 \\\\ 0 \\end{bmatrix} = \\begin{bmatrix} 3 \\\\ 0 \\end{bmatrix}$</button>
              <button class="m2-opt-btn" data-opt="D">$\\frac{15}{5}\\begin{bmatrix} 3 \\\\ 4 \\end{bmatrix} = \\begin{bmatrix} 9 \\\\ 12 \\end{bmatrix}$</button>
            </div>
            <div class="m2-fb" id="m2Fb2"></div>
          </div>

          <div class="m2-quiz-item" data-q="3" data-ans="D">
            <strong>Q3.</strong> Suppose $\\|\\mathbf{x}\\| = 3$ and $\\|\\mathbf{y}\\| = 4$. If $\\langle \\mathbf{x}, \\mathbf{y} \\rangle = -12$, what MUST be true?
            <div class="m2-opt-group">
              <button class="m2-opt-btn" data-opt="A">$\\mathbf{x}$ and $\\mathbf{y}$ are orthogonal.</button>
              <button class="m2-opt-btn" data-opt="B">$\\mathbf{x}$ and $\\mathbf{y}$ point in the exact same direction.</button>
              <button class="m2-opt-btn" data-opt="C">$\\|\\mathbf{x}+\\mathbf{y}\\| = 5$.</button>
              <button class="m2-opt-btn" data-opt="D">$\\mathbf{x}$ and $\\mathbf{y}$ point in exact opposite directions ($\\theta = 180^\\circ$).</button>
            </div>
            <div class="m2-fb" id="m2Fb3"></div>
          </div>
        </div>
      `;
    },

    mount(container, { LA, renderMath }) {
      // =====================================================================
      // LAB 1: 2D INTERACTIVE PROJECTION & CLOSEST-POINT CANVAS
      // =====================================================================
      const canvas = container.querySelector('#m2Canvas');
      const ctx = canvas.getContext('2d');
      const range = 5.5;

      let vx = [2, 3.5];
      let vy = [4, 1];
      let dragging = null;

      const snapEl = container.querySelector('#m2Snap');
      const tSlider = container.querySelector('#m2TSlider');
      const showCandEl = container.querySelector('#m2ShowCandidate');

      function getLayout() {
        const rect = canvas.getBoundingClientRect();
        const w = rect.width || 460;
        const h = rect.height || 370;
        const scale = Math.min(w, h) / (2 * range);
        return { w, h, cx: w / 2, cy: h / 2, scale };
      }

      function toScreen(pt, L) {
        return [L.cx + pt[0] * L.scale, L.cy - pt[1] * L.scale];
      }

      function toWorld(sx, sy, L) {
        let wx = (sx - L.cx) / L.scale;
        let wy = (L.cy - sy) / L.scale;
        if (snapEl.checked) {
          wx = Math.round(wx * 2) / 2;
          wy = Math.round(wy * 2) / 2;
        }
        return [Math.max(-5, Math.min(5, wx)), Math.max(-5, Math.min(5, wy))];
      }

      function drawArrow(x1, y1, x2, y2, color, width = 3, dashed = false) {
        const dx = x2 - x1, dy = y2 - y1;
        const len = Math.hypot(dx, dy);
        if (len < 2) return;
        const angle = Math.atan2(dy, dx);
        const head = 10;

        ctx.save();
        ctx.strokeStyle = color;
        ctx.fillStyle = color;
        ctx.lineWidth = width;
        if (dashed) ctx.setLineDash([6, 4]);

        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2 - head * 0.7 * Math.cos(angle), y2 - head * 0.7 * Math.sin(angle));
        ctx.stroke();

        ctx.setLineDash([]);
        ctx.beginPath();
        ctx.moveTo(x2, y2);
        ctx.lineTo(x2 - head * Math.cos(angle - Math.PI / 7), y2 - head * Math.sin(angle - Math.PI / 7));
        ctx.lineTo(x2 - head * Math.cos(angle + Math.PI / 7), y2 - head * Math.sin(angle + Math.PI / 7));
        ctx.closePath();
        ctx.fill();
        ctx.restore();
      }

      function drawLab1() {
        if (!canvas.isConnected) return;
        const dpr = window.devicePixelRatio || 1;
        const L = getLayout();
        canvas.width = Math.round(L.w * dpr);
        canvas.height = Math.round(L.h * dpr);
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.clearRect(0, 0, L.w, L.h);

        // Grid
        ctx.strokeStyle = '#F1F5F9';
        ctx.lineWidth = 1;
        for (let g = -5; g <= 5; g++) {
          const [gx] = toScreen([g, 0], L);
          const [, gy] = toScreen([0, g], L);
          ctx.beginPath(); ctx.moveTo(gx, 0); ctx.lineTo(gx, L.h); ctx.stroke();
          ctx.beginPath(); ctx.moveTo(0, gy); ctx.lineTo(L.w, gy); ctx.stroke();
        }

        // Axes
        ctx.strokeStyle = '#CBD5E1';
        ctx.lineWidth = 1.4;
        ctx.beginPath(); ctx.moveTo(0, L.cy); ctx.lineTo(L.w, L.cy); ctx.stroke();
        ctx.beginPath(); ctx.moveTo(L.cx, 0); ctx.lineTo(L.cx, L.h); ctx.stroke();

        const yy = LA.dot(vy, vy);
        const xy = LA.dot(vx, vy);
        const validY = yy > 1e-6;
        const tStar = validY ? xy / yy : 0;
        const p = validY ? LA.scale(tStar, vy) : [0, 0];
        const w = sub(vx, p);

        const [sx, sy] = toScreen(vx, L);
        const [yx, yyS] = toScreen(vy, L);
        const [px, py] = toScreen(p, L);

        // Infinite span{y} line
        if (validY) {
          const uy = LA.scale(1 / Math.sqrt(yy), vy);
          const pA = toScreen(LA.scale(-12, uy), L);
          const pB = toScreen(LA.scale(12, uy), L);
          ctx.save();
          ctx.strokeStyle = 'rgba(220, 38, 38, 0.25)';
          ctx.lineWidth = 2;
          ctx.setLineDash([7, 5]);
          ctx.beginPath();
          ctx.moveTo(pA[0], pA[1]);
          ctx.lineTo(pB[0], pB[1]);
          ctx.stroke();
          ctx.restore();
        }

        // Right-angle square marker at p between span{y} and w
        const nw = LA.norm(w);
        if (validY && nw > 0.25) {
          const uy = LA.scale(1 / Math.sqrt(yy), vy);
          const uw = LA.scale(1 / nw, w);
          const dirSign = LA.dot(p, vy) >= 0 ? -1 : 1;
          const sq = 0.26;
          const c1 = toScreen(LA.add(p, LA.scale(dirSign * sq, uy)), L);
          const c2 = toScreen(LA.add(LA.add(p, LA.scale(dirSign * sq, uy)), LA.scale(sq, uw)), L);
          const c3 = toScreen(LA.add(p, LA.scale(sq, uw)), L);
          ctx.save();
          ctx.strokeStyle = '#64748B';
          ctx.lineWidth = 1.6;
          ctx.beginPath();
          ctx.moveTo(c1[0], c1[1]);
          ctx.lineTo(c2[0], c2[1]);
          ctx.lineTo(c3[0], c3[1]);
          ctx.stroke();
          ctx.restore();
        }

        // Candidate test point t*y for Closest-Point Discovery
        const tCand = parseFloat(tSlider.value);
        if (validY && showCandEl.checked) {
          const cand = LA.scale(tCand, vy);
          const [cxS, cyS] = toScreen(cand, L);
          ctx.save();
          ctx.strokeStyle = '#7C3AED';
          ctx.lineWidth = 2;
          ctx.setLineDash([4, 4]);
          ctx.beginPath();
          ctx.moveTo(cxS, cyS);
          ctx.lineTo(sx, sy);
          ctx.stroke();
          ctx.restore();

          ctx.beginPath();
          ctx.arc(cxS, cyS, 5.5, 0, Math.PI * 2);
          ctx.fillStyle = '#7C3AED';
          ctx.fill();
          ctx.fillStyle = '#6D28D9';
          ctx.font = 'bold 12px sans-serif';
          ctx.fillText(`t·y`, cxS + 7, cyS + 14);
        }

        // Residual w (from p to x)
        if (nw > 0.05) {
          drawArrow(px, py, sx, sy, '#D97706', 2.6, true);
          ctx.fillStyle = '#B45309';
          ctx.font = 'bold 12px sans-serif';
          ctx.fillText('w = x - p', (px + sx) / 2 + 8, (py + sy) / 2);
        }

        // Vector y (red)
        drawArrow(L.cx, L.cy, yx, yyS, '#DC2626', 3.0);

        // Projection p (teal) drawn over y so it's clearly visible
        if (LA.norm(p) > 0.05) {
          drawArrow(L.cx, L.cy, px, py, '#0F766E', 3.6);
          ctx.fillStyle = '#0F766E';
          ctx.font = 'bold 13px sans-serif';
          ctx.fillText('p', px + 8, py - 8);
        }

        // Vector x (blue)
        drawArrow(L.cx, L.cy, sx, sy, '#1D4ED8', 3.2);

        // Draggable handles
        [[sx, sy, '#1D4ED8', 'x'], [yx, yyS, '#DC2626', 'y']].forEach(([hx, hy, col, lbl]) => {
          ctx.beginPath();
          ctx.arc(hx, hy, 7, 0, Math.PI * 2);
          ctx.fillStyle = '#FFFFFF';
          ctx.fill();
          ctx.lineWidth = 2.5;
          ctx.strokeStyle = col;
          ctx.stroke();
          ctx.fillStyle = col;
          ctx.font = 'bold 14px sans-serif';
          ctx.fillText(lbl, hx + 10, hy - 8);
        });

        // Update Live Stats
        if (!validY) {
          container.querySelector('#m2ProjVal').textContent = 'Undefined (y = 0)';
          container.querySelector('#m2ResVal').textContent = '—';
          container.querySelector('#m2OrthoChk').textContent = '—';
          container.querySelector('#m2PythChk').textContent = 'Drag y away from (0,0)';
          container.querySelector('#m2ClosestInfo').textContent = 'Set y ≠ 0 to explore projection.';
          return;
        }

        container.querySelector('#m2ProjVal').textContent =
          `(${p[0].toFixed(2)}, ${p[1].toFixed(2)})  [t* = ${tStar.toFixed(2)}]`;
        container.querySelector('#m2ResVal').textContent =
          `(${w[0].toFixed(2)}, ${w[1].toFixed(2)})  [‖w‖ = ${nw.toFixed(2)}]`;

        const wyDot = LA.dot(w, vy);
        container.querySelector('#m2OrthoChk').innerHTML =
          `<span style="color:#15803D;">${Math.abs(wyDot) < 1e-9 ? '0.00' : wyDot.toFixed(2)} ✓ (w ⊥ y)</span>`;

        const np2 = LA.dot(p, p), nw2 = LA.dot(w, w), nx2 = LA.dot(vx, vx);
        container.querySelector('#m2PythChk').textContent =
          `‖p‖² + ‖w‖² = ${np2.toFixed(2)} + ${nw2.toFixed(2)} = ${(np2 + nw2).toFixed(2)} = ‖x‖² (${nx2.toFixed(2)})`;

        const candErr = LA.norm(sub(vx, LA.scale(tCand, vy)));
        const isMin = Math.abs(tCand - tStar) < 0.04;
        container.querySelector('#m2ClosestInfo').innerHTML =
          `At <span style="color:#7C3AED;">t = ${tCand.toFixed(2)}</span>: ‖x − t·y‖ = <strong>${candErr.toFixed(3)}</strong> ` +
          `&nbsp;vs.&nbsp; Min distance ‖x − p‖ = <strong style="color:#0F766E;">${nw.toFixed(3)}</strong> ` +
          (isMin ? `<span style="color:#15803D;">(★ Optimal Projection Hit!)</span>` : `<span style="color:#64748B;">(Optimal t* = ${tStar.toFixed(2)})</span>`);
      }

      // Pointer handlers
      function getPointerPos(e) {
        const rect = canvas.getBoundingClientRect();
        const cx = e.touches ? e.touches[0].clientX : e.clientX;
        const cy = e.touches ? e.touches[0].clientY : e.clientY;
        return [cx - rect.left, cy - rect.top];
      }
      function onDown(e) {
        const L = getLayout();
        const [px, py] = getPointerPos(e);
        const [sx, sy] = toScreen(vx, L);
        const [yx, yyS] = toScreen(vy, L);
        const dx = Math.hypot(px - sx, py - sy);
        const dy = Math.hypot(px - yx, py - yyS);
        if (dx < 22 || dy < 22) {
          dragging = dx <= dy ? 'x' : 'y';
          e.preventDefault();
        }
      }
      function onMove(e) {
        if (!dragging || !canvas.isConnected) return;
        e.preventDefault();
        const L = getLayout();
        const [px, py] = getPointerPos(e);
        const pt = toWorld(px, py, L);
        if (dragging === 'x') vx = pt;
        else vy = pt;
        drawLab1();
      }
      function onUp() { dragging = null; }

      canvas.addEventListener('mousedown', onDown);
      window.addEventListener('mousemove', onMove);
      window.addEventListener('mouseup', onUp);
      canvas.addEventListener('touchstart', onDown, { passive: false });
      canvas.addEventListener('touchmove', onMove, { passive: false });
      canvas.addEventListener('touchend', onUp);

      const lab1Presets = {
        acute: [[2, 3.5], [4, 1]],
        obtuse: [[-2.5, 3], [4, 1]],
        ortho: [[-1, 3], [3, 1]],
        parallel: [[3, 1.5], [2, 1]]
      };
      container.querySelectorAll('[data-preset]').forEach(btn => {
        btn.addEventListener('click', () => {
          const [a, b] = lab1Presets[btn.dataset.preset];
          vx = [...a]; vy = [...b];
          drawLab1();
        });
      });

      tSlider.addEventListener('input', drawLab1);
      showCandEl.addEventListener('change', drawLab1);
      container.querySelector('#m2SnapOptimal').addEventListener('click', () => {
        const yy = LA.dot(vy, vy);
        if (yy > 1e-6) {
          const tStar = Math.max(-2.5, Math.min(2.5, LA.dot(vx, vy) / yy));
          tSlider.value = tStar.toFixed(2);
          drawLab1();
        }
      });
      window.addEventListener('resize', drawLab1);

      // =====================================================================
      // LAB 2: 3D PLOTLY PROJECTION & STEP-BY-STEP DECOMPOSITION
      // =====================================================================
      const x3Inputs = Array.from(container.querySelectorAll('.m2-3d-x'));
      const y3Inputs = Array.from(container.querySelectorAll('.m2-3d-y'));
      const plotlyDiv = container.querySelector('#m2Plotly3D');
      const calc3DDiv = container.querySelector('#m2Calc3D');

      function read3D(inputs) {
        return inputs.map(i => {
          const v = parseFloat(i.value);
          return Number.isFinite(v) ? v : 0;
        });
      }

      function updateLab2() {
        const x = read3D(x3Inputs);
        const y = read3D(y3Inputs);
        const xy = LA.dot(x, y);
        const yy = LA.dot(y, y);

        if (yy < 1e-10) {
          calc3DDiv.innerHTML = `<p style="color:#991B1B;"><strong>Note:</strong> $\\mathbf{y} = \\mathbf{0}$ cannot be projected onto. Please enter a nonzero vector $\\mathbf{y}$.</p>`;
          renderMath(calc3DDiv);
          return;
        }

        const t = xy / yy;
        const p = LA.scale(t, y);
        const w = sub(x, p);

        // Exact numerator arrays for clean LaTeX when integer inputs are used
        const pNum = y.map(yi => xy * yi);
        const wNum = x.map((xi, i) => xi * yy - xy * y[i]);
        const tTex = fmtFrac(xy, yy);

        calc3DDiv.innerHTML = `
          <p style="margin:0 0 6px;"><strong>Step 1 — Compute Inner Products:</strong></p>
          $$\\langle \\mathbf{x}, \\mathbf{y} \\rangle = ${Number(xy.toFixed(4))}, \\qquad \\langle \\mathbf{y}, \\mathbf{y} \\rangle = \\|\\mathbf{y}\\|^2 = ${Number(yy.toFixed(4))}$$
          <p style="margin:8px 0 6px;"><strong>Step 2 — Orthogonal Projection $\\mathbf{p} = \\operatorname{proj}_{\\mathbf{y}}\\mathbf{x}$:</strong></p>
          $$\\mathbf{p} = \\frac{\\langle \\mathbf{x}, \\mathbf{y} \\rangle}{\\langle \\mathbf{y}, \\mathbf{y} \\rangle}\\mathbf{y} = ${tTex}\\begin{bmatrix} ${y.join(' \\\\ ')} \\end{bmatrix} = \\begin{bmatrix} ${fmtVecExact(pNum, yy)} \\end{bmatrix}$$
          <p style="margin:8px 0 6px;"><strong>Step 3 — Orthogonal Component $\\mathbf{w} = \\mathbf{x} - \\mathbf{p}$:</strong></p>
          $$\\mathbf{w} = \\mathbf{x} - \\mathbf{p} = \\begin{bmatrix} ${fmtVecExact(wNum, yy)} \\end{bmatrix}$$
          <p style="margin:8px 0 0; font-size:0.88rem; color:#15803D;">
            <strong>Checks:</strong> $\\langle \\mathbf{w}, \\mathbf{y} \\rangle = ${Math.abs(LA.dot(w, y)) < 1e-9 ? '0' : LA.dot(w, y).toFixed(4)}$, &nbsp;
            $\\|\\mathbf{p}\\|^2 + \\|\\mathbf{w}\\|^2 = ${(LA.dot(p, p) + LA.dot(w, w)).toFixed(3)} = \\|\\mathbf{x}\\|^2 = ${LA.dot(x, x).toFixed(3)}$ ✓
          </p>
        `;
        renderMath(calc3DDiv);

        // Render 3D Plotly Scene
        if (typeof Plotly === 'undefined') {
          plotlyDiv.innerHTML = `<div style="padding:24px; color:var(--text-muted); text-align:center;">3D plot requires Plotly.js (check internet connection). Algebraic decomposition is shown on the right.</div>`;
          return;
        }

        const spanLo = Math.min(-0.4, t - 0.5, 0) - 0.3;
        const spanHi = Math.max(1.4, t + 0.5, 1) + 0.3;
        const lineA = LA.scale(spanLo, y);
        const lineB = LA.scale(spanHi, y);

        const traces = [
          // Span{y} dashed line
          {
            type: 'scatter3d', mode: 'lines', name: 'span{y}',
            x: [lineA[0], lineB[0]], y: [lineA[1], lineB[1]], z: [lineA[2], lineB[2]],
            line: { color: 'rgba(220,38,38,0.3)', width: 3, dash: 'dash' }
          },
          // Vector x
          {
            type: 'scatter3d', mode: 'lines+markers+text', name: 'x',
            x: [0, x[0]], y: [0, x[1]], z: [0, x[2]],
            text: ['', 'x'], textposition: 'top center',
            textfont: { color: '#1D4ED8', size: 13 },
            line: { color: '#1D4ED8', width: 6 },
            marker: { size: [3, 6], color: '#1D4ED8' }
          },
          // Vector y
          {
            type: 'scatter3d', mode: 'lines+markers+text', name: 'y',
            x: [0, y[0]], y: [0, y[1]], z: [0, y[2]],
            text: ['', 'y'], textposition: 'top center',
            textfont: { color: '#DC2626', size: 13 },
            line: { color: '#DC2626', width: 5 },
            marker: { size: [3, 5], color: '#DC2626' }
          },
          // Projection p
          {
            type: 'scatter3d', mode: 'lines+markers+text', name: 'p = proj_y(x)',
            x: [0, p[0]], y: [0, p[1]], z: [0, p[2]],
            text: ['', 'p'], textposition: 'bottom center',
            textfont: { color: '#0F766E', size: 13 },
            line: { color: '#0F766E', width: 7 },
            marker: { size: [3, 6], color: '#0F766E' }
          },
          // Orthogonal residual w = x - p (drawn from tip of p to tip of x)
          {
            type: 'scatter3d', mode: 'lines+text', name: 'w = x - p (⊥ y)',
            x: [p[0], x[0]], y: [p[1], x[1]], z: [p[2], x[2]],
            text: ['', ''],
            line: { color: '#D97706', width: 5, dash: 'dot' }
          }
        ];

        const layout = {
          margin: { l: 0, r: 0, t: 10, b: 0 },
          legend: { orientation: 'h', y: 0.99, x: 0.02, font: { size: 11 } },
          scene: {
            aspectmode: 'data',
            xaxis: { title: 'x₁', zerolinecolor: '#94A3B8' },
            yaxis: { title: 'x₂', zerolinecolor: '#94A3B8' },
            zaxis: { title: 'x₃', zerolinecolor: '#94A3B8' },
            camera: { eye: { x: 1.55, y: 1.45, z: 1.15 } }
          }
        };

        Plotly.react(plotlyDiv, traces, layout, { responsive: true, displayModeBar: false });
      }

      [...x3Inputs, ...y3Inputs].forEach(inp => inp.addEventListener('input', updateLab2));

      function set3DPreset(xv, yv) {
        x3Inputs.forEach((el, i) => el.value = xv[i]);
        y3Inputs.forEach((el, i) => el.value = yv[i]);
        updateLab2();
      }
      container.querySelector('#m2PresetEx1').addEventListener('click', () => set3DPreset([1, 2, 3], [1, 1, 1]));
      container.querySelector('#m2PresetEx2').addEventListener('click', () => set3DPreset([4, 1, 2], [1, 2, 2]));
      container.querySelector('#m2PresetEx3').addEventListener('click', () => set3DPreset([-2, 3, 1], [2, 1, -1]));

      // =====================================================================
      // SECTION D: CAUCHY-SCHWARZ & COSINE SIMILARITY TESTER (R^4)
      // =====================================================================
      const csUInputs = Array.from(container.querySelectorAll('.m2-cs-u'));
      const csVInputs = Array.from(container.querySelectorAll('.m2-cs-v'));
      const csReadout = container.querySelector('#m2CsReadout');
      const csBar = container.querySelector('#m2CsBar');

      function updateCS() {
        const u = read3D(csUInputs);
        const v = read3D(csVInputs);
        const dot = LA.dot(u, v);
        const nu = LA.norm(u);
        const nv = LA.norm(v);
        const prod = nu * nv;

        if (prod < 1e-10) {
          csReadout.innerHTML = `One of the vectors is $\\mathbf{0}$: $|\\langle \\mathbf{u},\\mathbf{v}\\rangle| = 0 = \\|\\mathbf{u}\\|\\|\\mathbf{v}\\|$ (equality holds trivially; angle undefined).`;
          csBar.style.width = '100%';
          renderMath(csReadout);
          return;
        }

        const cosT = Math.max(-1, Math.min(1, dot / prod));
        const deg = (Math.acos(cosT) * 180) / Math.PI;
        const pct = Math.min(100, (Math.abs(dot) / prod) * 100);
        const eq = Math.abs(Math.abs(dot) - prod) < 1e-6;

        csReadout.innerHTML = `
          $|\\langle \\mathbf{u}, \\mathbf{v} \\rangle| = \\mathbf{${Math.abs(dot).toFixed(3)}} ${eq ? '=' : '<'} \\|\\mathbf{u}\\|\\|\\mathbf{v}\\| = (${nu.toFixed(3)})(${nv.toFixed(3)}) = \\mathbf{${prod.toFixed(3)}}$
          ${eq ? ' <span style="color:#15803D; font-weight:700;">★ Equality! (Linearly Dependent)</span>' : ''}<br/>
          $\\cos\\theta = \\mathbf{${cosT.toFixed(4)}} \\quad \\Longrightarrow \\quad \\theta = \\mathbf{${deg.toFixed(1)}^\\circ}$
        `;
        csBar.style.width = `${pct.toFixed(1)}%`;
        renderMath(csReadout);
      }

      [...csUInputs, ...csVInputs].forEach(el => el.addEventListener('input', updateCS));
      container.querySelectorAll('[data-cs]').forEach(btn => {
        btn.addEventListener('click', () => {
          const mode = btn.dataset.cs;
          let u = [1, 2, -1, 2], v = [2, 4, -2, 4];
          if (mode === 'anti') { u = [1, -2, 0, 3]; v = [-2, 4, 0, -6]; }
          else if (mode === 'ortho') { u = [1, 2, -1, 0]; v = [2, -1, 0, 5]; }
          else if (mode === 'random') {
            const rnd = () => Math.floor(Math.random() * 9) - 4;
            u = [rnd(), rnd(), rnd(), rnd()];
            v = [rnd(), rnd(), rnd(), rnd()];
          }
          csUInputs.forEach((el, i) => el.value = u[i]);
          csVInputs.forEach((el, i) => el.value = v[i]);
          updateCS();
        });
      });

      // =====================================================================
      // SECTION E: QUICK CONCEPT CHECK
      // =====================================================================
      const explanations = {
        1: 'Set ⟨x, y⟩ = 2(4) + k(-2) + (-3)(2) = 8 - 2k - 6 = 2 - 2k = 0, which gives k = 1.',
        2: '⟨x, y⟩ = 15 and ‖y‖² = 3² + 4² = 25 (don’t forget to square ‖y‖!). Thus proj_y(x) = (15/25) y = (3/5) [3, 4]ᵀ.',
        3: 'Here cos θ = ⟨x, y⟩ / (‖x‖‖y‖) = -12 / (3·4) = -1, so θ = 180°. Equality holds in Cauchy–Schwarz with opposite directions!'
      };

      container.querySelectorAll('.m2-quiz-item').forEach(item => {
        const qNum = item.dataset.q;
        const ans = item.dataset.ans;
        const fb = container.querySelector(`#m2Fb${qNum}`);
        item.querySelectorAll('.m2-opt-btn').forEach(btn => {
          btn.addEventListener('click', () => {
            item.querySelectorAll('.m2-opt-btn').forEach(b => b.classList.remove('correct', 'wrong'));
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

      // Initial renders
      drawLab1();
      updateLab2();
      updateCS();
    }
  });
})();
