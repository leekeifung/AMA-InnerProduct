/**
 * AMA1751 Linear Algebra — Topic 6: Inner Product
 * Module 1: Inner Product and Norm in R^n (Slides P. 2–4)
 */
(function () {
  // Helper: format square root nicely (e.g., 16 -> "4", 55 -> "\sqrt{55}")
  function fmtSqrt(sumSq) {
    const rounded = Math.round(sumSq);
    if (Math.abs(sumSq - rounded) < 1e-8 && rounded >= 0) {
      const root = Math.round(Math.sqrt(rounded));
      if (root * root === rounded) return `${root}`;
      return `\\sqrt{${rounded}}`;
    }
    return Math.sqrt(Math.max(0, sumSq)).toFixed(3);
  }

  // Helper: format a signed number inside a sum/product term
  function fmtFactor(n) {
    const s = Number.isInteger(n) ? `${n}` : Number(n.toFixed(2)).toString();
    return n < 0 ? `(${s})` : s;
  }

  AMA1751.registerModule(1, {
    render(meta) {
      return `
        <style>
          .m1-canvas-wrap {
            position: relative;
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 10px;
            overflow: hidden;
          }
          .m1-canvas {
            display: block;
            width: 100%;
            height: 360px;
            cursor: grab;
            touch-action: none;
          }
          .m1-canvas:active { cursor: grabbing; }
          .m1-controls {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            margin-bottom: 12px;
          }
          .m1-chip {
            padding: 6px 11px;
            border-radius: 999px;
            border: 1px solid var(--border);
            background: #F8FAFC;
            font-size: 0.82rem;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.15s;
          }
          .m1-chip:hover { background: var(--polyu-red-soft); border-color: #E8B4C4; color: var(--polyu-red-dark); }
          .m1-stat-grid {
            display: grid;
            grid-template-columns: repeat(2, minmax(0, 1fr));
            gap: 10px;
          }
          .m1-stat {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 10px 12px;
          }
          .m1-stat-label { font-size: 0.76rem; color: var(--text-muted); font-weight: 600; text-transform: uppercase; letter-spacing: 0.04em; }
          .m1-stat-val { font-size: 1.05rem; font-weight: 700; margin-top: 3px; color: #0F172A; }
          .m1-badge {
            display: inline-block;
            padding: 2px 9px;
            border-radius: 999px;
            font-size: 0.76rem;
            font-weight: 700;
            margin-left: 6px;
          }
          .m1-badge.pos { background: #DCFCE7; color: #166534; }
          .m1-badge.neg { background: #FEE2E2; color: #991B1B; }
          .m1-badge.zero { background: #DBEAFE; color: #1E40AF; }
          .m1-vec-row {
            display: flex;
            align-items: center;
            gap: 8px;
            flex-wrap: wrap;
            margin-bottom: 10px;
          }
          .m1-vec-input {
            width: 64px;
            padding: 6px 8px;
            border: 1px solid #CBD5E1;
            border-radius: 6px;
            text-align: center;
            font-weight: 600;
            font-size: 0.92rem;
          }
          .m1-vec-input:focus { outline: 2px solid var(--polyu-red); border-color: var(--polyu-red); }
          .m1-math-scroll { overflow-x: auto; padding: 6px 2px; }
          .m1-prop-table {
            width: 100%;
            border-collapse: collapse;
            font-size: 0.9rem;
          }
          .m1-prop-table th, .m1-prop-table td {
            padding: 9px 10px;
            border-bottom: 1px solid var(--border);
            text-align: left;
          }
          .m1-prop-table th { background: #F8FAFC; color: #475569; font-size: 0.8rem; text-transform: uppercase; }
          .m1-quiz-item {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 14px 16px;
            margin-bottom: 12px;
          }
          .m1-opt-group { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 10px; }
          .m1-opt-btn {
            padding: 7px 13px;
            border: 1px solid #CBD5E1;
            border-radius: 8px;
            background: #FFFFFF;
            cursor: pointer;
            font-weight: 600;
            font-size: 0.88rem;
          }
          .m1-opt-btn:hover { background: #F1F5F9; }
          .m1-opt-btn.correct { background: #DCFCE7; border-color: #16A34A; color: #166534; }
          .m1-opt-btn.wrong { background: #FEE2E2; border-color: #DC2626; color: #991B1B; }
          .m1-fb { margin-top: 8px; font-size: 0.88rem; font-weight: 500; }
        </style>

        <!-- HEADER -->
        <div class="module-header">
          <span class="slide-ref">Part 1 · Lecture Notes ${meta.slides}</span>
          <h2>1. Inner Product and Norm in $\\mathbb{R}^n$</h2>
          <p>How we bring <strong>length</strong>, <strong>angle</strong>, and <strong>direction</strong> into $n$-dimensional space — the geometric foundation for orthogonality, projections, least squares, and the Spectral Theorem.</p>
        </div>

        <!-- SECTION A: PRACTICAL MOTIVATION & DEFINITIONS -->
        <div class="card">
          <h3>Why Do We Need an Inner Product &amp; Norm?</h3>
          <p>Up to now in AMA1751, you worked with vectors algebraically (linear combinations, spans, bases, solving $A\\mathbf{x}=\\mathbf{b}$). However, in engineering, data science, and physics, we also need <strong>geometry</strong>:</p>
          <ul>
            <li><strong>How long is a vector</strong> (or how large is a prediction error $\\mathbf{b} - A\\mathbf{x}$)? $\\longrightarrow$ Measured by the <strong>norm</strong> $\\|\\mathbf{x}\\|$.</li>
            <li><strong>Are two vectors pointing in a similar, opposite, or perpendicular direction?</strong> $\\longrightarrow$ Measured by the <strong>inner product</strong> $\\langle \\mathbf{x}, \\mathbf{y} \\rangle$.</li>
            <li><strong>How do we extract pure direction without magnitude?</strong> $\\longrightarrow$ Scale to a <strong>unit vector</strong> $\\mathbf{u} = \\frac{1}{\\|\\mathbf{x}\\|}\\mathbf{x}$.</li>
          </ul>

          <div class="grid-2">
            <div class="callout definition">
              <strong>Definition — Inner Product (Scalar Product):</strong><br/>
              Given $\\mathbf{x} = \\begin{bmatrix} x_1 \\\\ \\vdots \\\\ x_n \\end{bmatrix}$ and $\\mathbf{y} = \\begin{bmatrix} y_1 \\\\ \\vdots \\\\ y_n \\end{bmatrix}$ in $\\mathbb{R}^n$, their inner product is the scalar:
              $$\\langle \\mathbf{x}, \\mathbf{y} \\rangle = \\mathbf{y}^T \\mathbf{x} = \\begin{bmatrix} y_1 & y_2 & \\cdots & y_n \\end{bmatrix}\\begin{bmatrix} x_1 \\\\ x_2 \\\\ \\vdots \\\\ x_n \\end{bmatrix} = x_1 y_1 + x_2 y_2 + \\cdots + x_n y_n.$$
            </div>

            <div class="callout definition">
              <strong>Definition — Norm &amp; Unit Vector:</strong><br/>
              1. The <strong>norm</strong> (Euclidean length) of $\\mathbf{x} \\in \\mathbb{R}^n$ is:
              $$\\|\\mathbf{x}\\| = \\sqrt{\\langle \\mathbf{x}, \\mathbf{x} \\rangle} = \\sqrt{x_1^2 + x_2^2 + \\cdots + x_n^2}.$$
              2. A vector $\\mathbf{u}$ is called a <strong>unit vector</strong> if $\\|\\mathbf{u}\\| = 1$.<br/>
              <span style="font-size:0.88rem; color:var(--text-muted);">For any nonzero $\\mathbf{x}$, $\\mathbf{z} = \\frac{1}{\\|\\mathbf{x}\\|}\\mathbf{x}$ is always a unit vector in the same direction.</span>
            </div>
          </div>
        </div>

        <!-- SECTION B: INTERACTIVE LAB 1 — 2D GEOMETRIC SANDBOX -->
        <div class="card">
          <h3>Interactive Lab 1: Geometric Intuition in $\\mathbb{R}^2$ (Drag the Vector Tips!)</h3>
          <p style="margin-top:0; color:var(--text-muted);">Drag the <strong style="color:#1D4ED8;">blue tip ($\\mathbf{x}$)</strong> or <strong style="color:#DC2626;">red tip ($\\mathbf{y}$)</strong> directly on the plane to see how the inner product sign, angle, scaling, and the Triangle Inequality behave in real time.</p>

          <div class="m1-controls">
            <button class="m1-chip" data-preset="acute">Acute ($\\langle \\mathbf{x},\\mathbf{y}\\rangle > 0$)</button>
            <button class="m1-chip" data-preset="obtuse">Obtuse ($\\langle \\mathbf{x},\\mathbf{y}\\rangle < 0$)</button>
            <button class="m1-chip" data-preset="ortho">Perpendicular ($= 0$)</button>
            <button class="m1-chip" data-preset="parallel">Parallel (Triangle Equality)</button>
            <label class="m1-chip" style="display:inline-flex; align-items:center; gap:5px;">
              <input type="checkbox" id="m1ShowTri" checked /> Show $\\mathbf{x}+\\mathbf{y}$ Triangle
            </label>
            <label class="m1-chip" style="display:inline-flex; align-items:center; gap:5px;">
              <input type="checkbox" id="m1ShowUnit" /> Show Unit Circle &amp; $\\hat{\\mathbf{x}}$
            </label>
            <label class="m1-chip" style="display:inline-flex; align-items:center; gap:5px;">
              <input type="checkbox" id="m1Snap" checked /> Snap to 0.5
            </label>
          </div>

          <div class="grid-2" style="align-items:start;">
            <div class="m1-canvas-wrap">
              <canvas id="m1Canvas" class="m1-canvas"></canvas>
            </div>

            <div>
              <div class="m1-stat-grid">
                <div class="m1-stat">
                  <div class="m1-stat-label">Vectors</div>
                  <div class="m1-stat-val" id="m1VecText" style="font-size:0.92rem;"></div>
                </div>
                <div class="m1-stat">
                  <div class="m1-stat-label">Inner Product $\\langle \\mathbf{x}, \\mathbf{y} \\rangle$</div>
                  <div class="m1-stat-val" id="m1DotText"></div>
                </div>
                <div class="m1-stat">
                  <div class="m1-stat-label">Norms $\\|\\mathbf{x}\\|, \\|\\mathbf{y}\\|$</div>
                  <div class="m1-stat-val" id="m1NormText" style="font-size:0.94rem;"></div>
                </div>
                <div class="m1-stat">
                  <div class="m1-stat-label">Angle $\\theta$</div>
                  <div class="m1-stat-val" id="m1AngleText"></div>
                </div>
              </div>

              <!-- Triangle Inequality Live Meter -->
              <div class="callout tip" style="margin-top:12px;">
                <strong>Triangle Inequality Check (Slide P. 4):</strong>
                <div id="m1TriText" style="margin-top:4px; font-weight:600;"></div>
                <div style="margin-top:8px; background:#E2E8F0; height:10px; border-radius:999px; overflow:hidden;">
                  <div id="m1TriBar" style="height:100%; background:var(--accent-teal); width:70%; transition:width 0.1s;"></div>
                </div>
                <div style="font-size:0.78rem; color:var(--text-muted); margin-top:4px;">
                  Bar shows ratio $\\|\\mathbf{x}+\\mathbf{y}\\| / (\\|\\mathbf{x}\\| + \\|\\mathbf{y}\\|) \\le 100\\%$. Try making $\\mathbf{x}$ and $\\mathbf{y}$ point in the exact same direction!
                </div>
              </div>

              <!-- Scalar Scaling Slider (Slide P. 4: ||cx|| = |c| ||x||) -->
              <div class="m1-stat" style="margin-top:12px;">
                <div style="display:flex; justify-content:space-between; align-items:center;">
                  <span class="m1-stat-label">Scalar Multiplier $c$ (test negative $c$!)</span>
                  <strong id="m1ScalarVal" style="color:var(--polyu-red);">c = 1.0</strong>
                </div>
                <input type="range" id="m1ScalarSlider" min="-2.5" max="2.5" step="0.5" value="1" style="width:100%; margin:8px 0;" />
                <div id="m1ScalarInfo" style="font-size:0.86rem; color:#334155;"></div>
              </div>
            </div>
          </div>
        </div>

        <!-- SECTION C: INTERACTIVE LAB 2 — R^n CALCULATOR & SLIDE P. 3 WORKED EXAMPLE -->
        <div class="card">
          <h3>Interactive Lab 2: Worked Example (Slide P. 3) &amp; $\\mathbb{R}^n$ Calculator</h3>
          <p style="margin-top:0; color:var(--text-muted);">Pre-loaded with the lecture example in $\\mathbb{R}^5$. Edit any coordinate or change dimension $n$ to see the complete step-by-step inner product, norm, and unit vector calculation.</p>

          <div class="m1-controls">
            <label style="font-size:0.88rem; font-weight:600; display:inline-flex; align-items:center; gap:6px;">
              Dimension $n$:
              <select id="m1DimSelect" style="padding:5px 10px; border-radius:6px; border:1px solid #CBD5E1; font-weight:600;">
                <option value="2">ℝ² (n = 2)</option>
                <option value="3">ℝ³ (n = 3)</option>
                <option value="4">ℝ⁴ (n = 4)</option>
                <option value="5" selected>ℝ⁵ (n = 5, Slide P. 3)</option>
              </select>
            </label>
            <button class="btn" id="m1LoadSlide3">Reset to Slide P. 3 Example</button>
            <button class="btn" id="m1RandomVecs">Random Integer Vectors</button>
          </div>

          <div class="m1-vec-row">
            <strong style="width:52px; color:#1D4ED8;">$\\mathbf{x}^T =$</strong>
            <div id="m1XInputs" style="display:flex; gap:6px; flex-wrap:wrap;"></div>
          </div>
          <div class="m1-vec-row">
            <strong style="width:52px; color:#DC2626;">$\\mathbf{y}^T =$</strong>
            <div id="m1YInputs" style="display:flex; gap:6px; flex-wrap:wrap;"></div>
          </div>

          <div id="m1CalcOutput" class="m1-math-scroll" style="background:#F8FAFC; border:1px solid var(--border); border-radius:10px; padding:14px 18px; margin-top:12px;"></div>

          <details class="code-details">
            <summary>Python (SymPy) Companion Code — Matches Lecture Slide P. 3</summary>
            <pre><code>import sympy as sp

x = sp.Matrix([[1, 2, 3, 4, 5]]).T
y = sp.Matrix([[3, -1, 2, 1, -1]]).T

# Inner product of x and y: y^T x
print("Inner product:", y.T @ x)   # Matrix([[6]])

# Norms of x and y
print("||x|| =", x.norm())         # sqrt(55)
print("||y|| =", y.norm())         # 4

# Normalize x into a unit vector z
z = x / x.norm()
print("Unit vector z =", z)
print("||z|| =", z.norm())         # 1</code></pre>
          </details>
        </div>

        <!-- SECTION D: PROPERTIES OF INNER PRODUCT & NORM (SLIDE P. 4) -->
        <div class="card">
          <h3>Properties of Inner Product &amp; Norm — Live Verifier (Slide P. 4)</h3>
          <p style="margin-top:0; color:var(--text-muted);">Verify every algebraic property on real vectors in $\\mathbb{R}^3$. Click the button below to test fresh random vectors and negative/positive scalars.</p>

          <div class="m1-controls" style="align-items:center;">
            <button class="btn btn-primary" id="m1NewPropTest">Generate New Test Vectors $(\\mathbf{x}, \\mathbf{y}, \\mathbf{w}, c)$</button>
            <span id="m1PropVectorsLabel" style="font-size:0.86rem; color:#334155;"></span>
          </div>

          <div class="grid-2">
            <div>
              <h4 style="margin:6px 0 8px; color:var(--accent-blue);">Properties of Inner Product</h4>
              <table class="m1-prop-table" id="m1IpPropTable"></table>
            </div>
            <div>
              <h4 style="margin:6px 0 8px; color:var(--accent-teal);">Properties of Norm</h4>
              <table class="m1-prop-table" id="m1NormPropTable"></table>
            </div>
          </div>
        </div>

        <!-- SECTION E: QUICK CONCEPT CHECK -->
        <div class="card">
          <h3>Quick Concept Check (Before Moving to Orthogonality)</h3>

          <div class="m1-quiz-item" data-q="1" data-ans="B">
            <strong>Q1.</strong> Let $\\mathbf{x} = \\begin{bmatrix} 2 \\\\ -3 \\\\ 1 \\end{bmatrix}$ and $\\mathbf{y} = \\begin{bmatrix} 1 \\\\ 2 \\\\ -1 \\end{bmatrix}$. What is $\\langle \\mathbf{x}, \\mathbf{y} \\rangle$, and is the angle between them acute or obtuse?
            <div class="m1-opt-group">
              <button class="m1-opt-btn" data-opt="A">$5$ (acute)</button>
              <button class="m1-opt-btn" data-opt="B">$-5$ (obtuse)</button>
              <button class="m1-opt-btn" data-opt="C">$-3$ (obtuse)</button>
              <button class="m1-opt-btn" data-opt="D">$0$ (right angle)</button>
            </div>
            <div class="m1-fb" id="m1Fb1"></div>
          </div>

          <div class="m1-quiz-item" data-q="2" data-ans="C">
            <strong>Q2.</strong> Suppose $\\|\\mathbf{x}\\| = 5$. What is the value of $\\|-3\\mathbf{x}\\|$?
            <div class="m1-opt-group">
              <button class="m1-opt-btn" data-opt="A">$-15$</button>
              <button class="m1-opt-btn" data-opt="B">$8$</button>
              <button class="m1-opt-btn" data-opt="C">$15$</button>
              <button class="m1-opt-btn" data-opt="D">$45$</button>
            </div>
            <div class="m1-fb" id="m1Fb2"></div>
          </div>

          <div class="m1-quiz-item" data-q="3" data-ans="A">
            <strong>Q3.</strong> Which vector is a <strong>unit vector</strong> in the same direction as $\\mathbf{v} = \\begin{bmatrix} 3 \\\\ -4 \\\\ 0 \\end{bmatrix}$?
            <div class="m1-opt-group">
              <button class="m1-opt-btn" data-opt="A">$\\frac{1}{5}\\begin{bmatrix} 3 \\\\ -4 \\\\ 0 \\end{bmatrix}$</button>
              <button class="m1-opt-btn" data-opt="B">$\\frac{1}{25}\\begin{bmatrix} 3 \\\\ -4 \\\\ 0 \\end{bmatrix}$</button>
              <button class="m1-opt-btn" data-opt="C">$\\frac{1}{7}\\begin{bmatrix} 3 \\\\ -4 \\\\ 0 \\end{bmatrix}$</button>
              <button class="m1-opt-btn" data-opt="D">$\\begin{bmatrix} 1 \\\\ -1 \\\\ 0 \\end{bmatrix}$</button>
            </div>
            <div class="m1-fb" id="m1Fb3"></div>
          </div>
        </div>
      `;
    },

    mount(container, { LA, renderMath }) {
      // =====================================================================
      // LAB 1: INTERACTIVE 2D GEOMETRIC CANVAS
      // =====================================================================
      const canvas = container.querySelector('#m1Canvas');
      const ctx = canvas.getContext('2d');
      const range = 5.5; // Coordinate bounds [-5.5, 5.5]

      let vx = [3, 1];
      let vy = [1.5, 3];
      let scalarC = 1.0;
      let dragging = null; // 'x' or 'y'

      const showTriEl = container.querySelector('#m1ShowTri');
      const showUnitEl = container.querySelector('#m1ShowUnit');
      const snapEl = container.querySelector('#m1Snap');
      const sliderEl = container.querySelector('#m1ScalarSlider');

      function getLayout() {
        const rect = canvas.getBoundingClientRect();
        const w = rect.width || 460;
        const h = rect.height || 360;
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
        wx = Math.max(-5, Math.min(5, wx));
        wy = Math.max(-5, Math.min(5, wy));
        return [wx, wy];
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

      function drawCanvas() {
        const dpr = window.devicePixelRatio || 1;
        const L = getLayout();
        canvas.width = Math.round(L.w * dpr);
        canvas.height = Math.round(L.h * dpr);
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

        ctx.clearRect(0, 0, L.w, L.h);

        // Grid lines
        ctx.strokeStyle = '#F1F5F9';
        ctx.lineWidth = 1;
        for (let g = -5; g <= 5; g++) {
          const [gx] = toScreen([g, 0], L);
          const [, gy] = toScreen([0, g], L);
          ctx.beginPath(); ctx.moveTo(gx, 0); ctx.lineTo(gx, L.h); ctx.stroke();
          ctx.beginPath(); ctx.moveTo(0, gy); ctx.lineTo(L.w, gy); ctx.stroke();
        }

        // Axes
        ctx.strokeStyle = '#94A3B8';
        ctx.lineWidth = 1.5;
        ctx.beginPath(); ctx.moveTo(0, L.cy); ctx.lineTo(L.w, L.cy); ctx.stroke();
        ctx.beginPath(); ctx.moveTo(L.cx, 0); ctx.lineTo(L.cx, L.h); ctx.stroke();

        // Axis ticks
        ctx.fillStyle = '#64748B';
        ctx.font = '11px sans-serif';
        for (let g = -5; g <= 5; g++) {
          if (g === 0) continue;
          const [gx] = toScreen([g, 0], L);
          const [, gy] = toScreen([0, g], L);
          ctx.fillText(g, gx - 3, L.cy + 13);
          ctx.fillText(g, L.cx + 5, gy + 4);
        }

        // Optional Unit Circle & normalized x_hat
        if (showUnitEl.checked) {
          ctx.save();
          ctx.strokeStyle = '#8B1E3F';
          ctx.lineWidth = 1.5;
          ctx.setLineDash([4, 4]);
          ctx.beginPath();
          ctx.arc(L.cx, L.cy, L.scale, 0, Math.PI * 2);
          ctx.stroke();
          ctx.restore();

          if (LA.norm(vx) > 1e-6) {
            const u = LA.normalize(vx);
            const [ux, uy] = toScreen(u, L);
            drawArrow(L.cx, L.cy, ux, uy, '#8B1E3F', 2.5);
            ctx.fillStyle = '#8B1E3F';
            ctx.font = 'bold 12px sans-serif';
            ctx.fillText('x̂ (||x̂||=1)', ux + 6, uy - 4);
          }
        }

        // Angle arc between x and y
        const nx = LA.norm(vx), ny = LA.norm(vy);
        if (nx > 0.1 && ny > 0.1) {
          const a1 = Math.atan2(-vx[1], vx[0]);
          const a2 = Math.atan2(-vy[1], vy[0]);
          let diff = a2 - a1;
          while (diff > Math.PI) diff -= 2 * Math.PI;
          while (diff < -Math.PI) diff += 2 * Math.PI;
          ctx.save();
          ctx.beginPath();
          ctx.moveTo(L.cx, L.cy);
          ctx.arc(L.cx, L.cy, 26, a1, a1 + diff, diff < 0);
          ctx.closePath();
          const dot = LA.dot(vx, vy);
          ctx.fillStyle = Math.abs(dot) < 1e-6 ? 'rgba(29,78,216,0.18)' : (dot > 0 ? 'rgba(22,163,74,0.16)' : 'rgba(220,38,38,0.16)');
          ctx.fill();
          ctx.restore();
        }

        const [sx, sy] = toScreen(vx, L);
        const [tx, ty] = toScreen(vy, L);

        // Optional Triangle Inequality overlay (x + y)
        if (showTriEl.checked) {
          const sum = LA.add(vx, vy);
          const [sumX, sumY] = toScreen(sum, L);
          drawArrow(sx, sy, sumX, sumY, '#15803D', 2, true); // shifted y
          drawArrow(L.cx, L.cy, sumX, sumY, '#0F766E', 2.6); // x + y
          ctx.fillStyle = '#0F766E';
          ctx.font = 'bold 12px sans-serif';
          ctx.fillText('x + y', sumX + 6, sumY - 6);
        }

        // Scaled vector c*x if c != 1
        if (Math.abs(scalarC - 1) > 1e-6) {
          const cxv = LA.scale(scalarC, vx);
          const [cxS, cyS] = toScreen(cxv, L);
          drawArrow(L.cx, L.cy, cxS, cyS, '#C59B27', 2.6, true);
          ctx.fillStyle = '#92400E';
          ctx.font = 'bold 12px sans-serif';
          ctx.fillText(`${scalarC}x`, cxS + 6, cyS + 14);
        }

        // Primary vectors x (blue) and y (red)
        drawArrow(L.cx, L.cy, sx, sy, '#1D4ED8', 3.2);
        drawArrow(L.cx, L.cy, tx, ty, '#DC2626', 3.2);

        // Draggable handles
        [[sx, sy, '#1D4ED8', 'x'], [tx, ty, '#DC2626', 'y']].forEach(([px, py, col, lbl]) => {
          ctx.beginPath();
          ctx.arc(px, py, 7, 0, Math.PI * 2);
          ctx.fillStyle = '#FFFFFF';
          ctx.fill();
          ctx.lineWidth = 2.5;
          ctx.strokeStyle = col;
          ctx.stroke();

          ctx.fillStyle = col;
          ctx.font = 'bold 14px sans-serif';
          ctx.fillText(lbl, px + 10, py - 8);
        });

        updateLab1Stats();
      }

      function updateLab1Stats() {
        const dot = LA.dot(vx, vy);
        const nx = LA.norm(vx);
        const ny = LA.norm(vy);
        const sumNorm = LA.norm(LA.add(vx, vy));

        container.querySelector('#m1VecText').innerHTML =
          `<span style="color:#1D4ED8;">x = (${vx[0].toFixed(1)}, ${vx[1].toFixed(1)})</span><br/>` +
          `<span style="color:#DC2626;">y = (${vy[0].toFixed(1)}, ${vy[1].toFixed(1)})</span>`;

        let badge = '<span class="m1-badge zero">Orthogonal (= 0)</span>';
        if (dot > 1e-6) badge = '<span class="m1-badge pos">Acute (> 0)</span>';
        else if (dot < -1e-6) badge = '<span class="m1-badge neg">Obtuse (< 0)</span>';
        container.querySelector('#m1DotText').innerHTML = `${dot.toFixed(2)} ${badge}`;

        container.querySelector('#m1NormText').innerHTML =
          `‖x‖ = ${nx.toFixed(2)} &nbsp;|&nbsp; ‖y‖ = ${ny.toFixed(2)}`;

        if (nx > 1e-6 && ny > 1e-6) {
          const cosT = Math.max(-1, Math.min(1, dot / (nx * ny)));
          const deg = (Math.acos(cosT) * 180) / Math.PI;
          container.querySelector('#m1AngleText').innerHTML =
            `${deg.toFixed(1)}° <span style="font-size:0.78rem; color:var(--text-muted); font-weight:500;">(cos θ = ${cosT.toFixed(2)})</span>`;
        } else {
          container.querySelector('#m1AngleText').textContent = 'Undefined (zero vector)';
        }

        const sumIndiv = nx + ny;
        const ratio = sumIndiv > 1e-6 ? Math.min(100, (sumNorm / sumIndiv) * 100) : 0;
        container.querySelector('#m1TriText').textContent =
          `‖x + y‖ = ${sumNorm.toFixed(2)}  ≤  ‖x‖ + ‖y‖ = ${nx.toFixed(2)} + ${ny.toFixed(2)} = ${sumIndiv.toFixed(2)}`;
        container.querySelector('#m1TriBar').style.width = `${ratio}%`;

        // Scalar slider check
        const scaledNorm = LA.norm(LA.scale(scalarC, vx));
        const scaledDot = LA.dot(LA.scale(scalarC, vx), vy);
        container.querySelector('#m1ScalarVal').textContent = `c = ${scalarC.toFixed(1)}`;
        container.querySelector('#m1ScalarInfo').innerHTML =
          `<strong>Norm:</strong> ‖(${scalarC.toFixed(1)})x‖ = <strong>${scaledNorm.toFixed(2)}</strong> = |${scalarC.toFixed(1)}| · ${nx.toFixed(2)}<br/>` +
          `<strong>Inner product:</strong> ⟨(${scalarC.toFixed(1)})x, y⟩ = <strong>${scaledDot.toFixed(2)}</strong> = (${scalarC.toFixed(1)}) · (${dot.toFixed(2)})`;
      }

      // Pointer Events (Mouse + Touch)
      function getPointerPos(e) {
        const rect = canvas.getBoundingClientRect();
        const clientX = e.touches ? e.touches[0].clientX : e.clientX;
        const clientY = e.touches ? e.touches[0].clientY : e.clientY;
        return [clientX - rect.left, clientY - rect.top];
      }

      function onDown(e) {
        const L = getLayout();
        const [px, py] = getPointerPos(e);
        const [sx, sy] = toScreen(vx, L);
        const [tx, ty] = toScreen(vy, L);
        const dx = Math.hypot(px - sx, py - sy);
        const dy = Math.hypot(px - tx, py - ty);
        if (dx < 22 || dy < 22) {
          dragging = dx <= dy ? 'x' : 'y';
          e.preventDefault();
        }
      }

      function onMove(e) {
        if (!dragging) return;
        e.preventDefault();
        const L = getLayout();
        const [px, py] = getPointerPos(e);
        const worldPt = toWorld(px, py, L);
        if (dragging === 'x') vx = worldPt;
        else vy = worldPt;
        drawCanvas();
      }

      function onUp() { dragging = null; }

      canvas.addEventListener('mousedown', onDown);
      window.addEventListener('mousemove', onMove);
      window.addEventListener('mouseup', onUp);
      canvas.addEventListener('touchstart', onDown, { passive: false });
      canvas.addEventListener('touchmove', onMove, { passive: false });
      canvas.addEventListener('touchend', onUp);

      // Presets & Toggles
      const presets = {
        acute: [[3.5, 0.5], [1.5, 3]],
        obtuse: [[3.5, 0], [-1.5, 3]],
        ortho: [[3, 1], [-1, 3]],
        parallel: [[2, 1], [3, 1.5]]
      };
      container.querySelectorAll('[data-preset]').forEach(btn => {
        btn.addEventListener('click', () => {
          const [px, py] = presets[btn.dataset.preset];
          vx = [...px]; vy = [...py];
          drawCanvas();
        });
      });
      showTriEl.addEventListener('change', drawCanvas);
      showUnitEl.addEventListener('change', drawCanvas);
      sliderEl.addEventListener('input', () => {
        scalarC = parseFloat(sliderEl.value);
        drawCanvas();
      });
      window.addEventListener('resize', drawCanvas);

      // =====================================================================
      // LAB 2: R^n CALCULATOR & SLIDE P. 3 WORKED EXAMPLE
      // =====================================================================
      const dimSelect = container.querySelector('#m1DimSelect');
      const xInputsDiv = container.querySelector('#m1XInputs');
      const yInputsDiv = container.querySelector('#m1YInputs');
      const calcOut = container.querySelector('#m1CalcOutput');

      let calcX = [1, 2, 3, 4, 5];
      let calcY = [3, -1, 2, 1, -1];

      function buildVectorInputs() {
        const n = parseInt(dimSelect.value, 10);
        while (calcX.length < n) { calcX.push(1); calcY.push(0); }
        calcX = calcX.slice(0, n);
        calcY = calcY.slice(0, n);

        xInputsDiv.innerHTML = calcX.map((v, i) =>
          `<input type="number" step="any" class="m1-vec-input" data-vec="x" data-idx="${i}" value="${v}" />`
        ).join('');
        yInputsDiv.innerHTML = calcY.map((v, i) =>
          `<input type="number" step="any" class="m1-vec-input" data-vec="y" data-idx="${i}" value="${v}" />`
        ).join('');

        container.querySelectorAll('.m1-vec-input').forEach(inp => {
          inp.addEventListener('input', () => {
            const idx = parseInt(inp.dataset.idx, 10);
            const val = parseFloat(inp.value);
            const safeVal = Number.isFinite(val) ? val : 0;
            if (inp.dataset.vec === 'x') calcX[idx] = safeVal;
            else calcY[idx] = safeVal;
            renderLab2Calc();
          });
        });

        renderLab2Calc();
      }

      function renderLab2Calc() {
        const dot = LA.dot(calcX, calcY);
        const sx2 = LA.dot(calcX, calcX);
        const sy2 = LA.dot(calcY, calcY);
        const nx = Math.sqrt(sx2);
        const ny = Math.sqrt(sy2);

        const dotTerms = calcX.map((xi, i) => `${fmtFactor(xi)}\\cdot ${fmtFactor(calcY[i])}`).join(' + ');
        const xSqTerms = calcX.map(xi => `${fmtFactor(xi)}^2`).join(' + ');
        const ySqTerms = calcY.map(yi => `${fmtFactor(yi)}^2`).join(' + ');

        const nxTex = fmtSqrt(sx2);
        const nyTex = fmtSqrt(sy2);

        let unitBlock = '';
        if (nx > 1e-10) {
          const z = LA.normalize(calcX);
          const zDec = z.map(v => v.toFixed(4)).join(' \\\\ ');
          unitBlock = `
            <p style="margin:10px 0 0;"><strong>Unit Vector in the direction of $\\mathbf{x}$:</strong></p>
            $$\\mathbf{z} = \\frac{1}{\\|\\mathbf{x}\\|}\\mathbf{x} = \\frac{1}{${nxTex}} \\begin{bmatrix} ${calcX.join(' \\\\ ')} \\end{bmatrix} \\approx \\begin{bmatrix} ${zDec} \\end{bmatrix}, \\qquad \\|\\mathbf{z}\\| = ${LA.norm(z).toFixed(4)} = 1.$$
          `;
        } else {
          unitBlock = `<p style="color:#991B1B; margin-top:10px;"><strong>Note:</strong> $\\mathbf{x} = \\mathbf{0}$ has norm $0$ and cannot be normalized into a unit vector.</p>`;
        }

        calcOut.innerHTML = `
          $$\\langle \\mathbf{x}, \\mathbf{y} \\rangle = \\mathbf{y}^T \\mathbf{x} = ${dotTerms} = ${Number(dot.toFixed(4))}$$
          $$\\|\\mathbf{x}\\| = \\sqrt{${xSqTerms}} = ${nxTex} ${nxTex.includes('sqrt') ? `\\approx ${nx.toFixed(4)}` : ''}$$
          $$\\|\\mathbf{y}\\| = \\sqrt{${ySqTerms}} = ${nyTex} ${nyTex.includes('sqrt') ? `\\approx ${ny.toFixed(4)}` : ''}$$
          ${unitBlock}
        `;
        renderMath(calcOut);
      }

      dimSelect.addEventListener('change', buildVectorInputs);
      container.querySelector('#m1LoadSlide3').addEventListener('click', () => {
        dimSelect.value = '5';
        calcX = [1, 2, 3, 4, 5];
        calcY = [3, -1, 2, 1, -1];
        buildVectorInputs();
      });
      container.querySelector('#m1RandomVecs').addEventListener('click', () => {
        const n = parseInt(dimSelect.value, 10);
        const rnd = () => Math.floor(Math.random() * 9) - 4; // -4..4
        calcX = Array.from({ length: n }, rnd);
        calcY = Array.from({ length: n }, rnd);
        buildVectorInputs();
      });

      // =====================================================================
      // SECTION D: PROPERTIES VERIFIER (SLIDE P. 4)
      // =====================================================================
      function runPropVerifier() {
        const rnd = () => Math.floor(Math.random() * 7) - 3;
        const x = [rnd() || 2, rnd(), rnd() || -1];
        const y = [rnd(), rnd() || 1, rnd()];
        const w = [rnd() || 1, rnd(), rnd()];
        const scalars = [-3, -2, -1.5, 2, 3];
        const c = scalars[Math.floor(Math.random() * scalars.length)];

        container.querySelector('#m1PropVectorsLabel').innerHTML =
          `Testing $\\mathbf{x}=(${x.join(',')})^T$, $\\mathbf{y}=(${y.join(',')})^T$, $\\mathbf{w}=(${w.join(',')})^T$, $c=${c}$`;

        const ipRows = [
          ['1. $\\langle \\mathbf{x},\\mathbf{y}\\rangle = \\langle \\mathbf{y},\\mathbf{x}\\rangle$', LA.dot(x, y), LA.dot(y, x), '='],
          ['2. $\\langle \\mathbf{x}+\\mathbf{y},\\mathbf{w}\\rangle = \\langle \\mathbf{x},\\mathbf{w}\\rangle + \\langle \\mathbf{y},\\mathbf{w}\\rangle$', LA.dot(LA.add(x, y), w), LA.dot(x, w) + LA.dot(y, w), '='],
          ['3. $\\langle c\\mathbf{x},\\mathbf{y}\\rangle = c\\langle \\mathbf{x},\\mathbf{y}\\rangle$', LA.dot(LA.scale(c, x), y), c * LA.dot(x, y), '='],
          ['4. $\\langle \\mathbf{x},\\mathbf{x}\\rangle \\ge 0$', LA.dot(x, x), 0, '≥']
        ];

        const normRows = [
          ['1. $\\|c\\mathbf{x}\\| = |c|\\cdot\\|\\mathbf{x}\\|$', LA.norm(LA.scale(c, x)), Math.abs(c) * LA.norm(x), '='],
          ['2. $\\|\\mathbf{x}+\\mathbf{y}\\| \\le \\|\\mathbf{x}\\| + \\|\\mathbf{y}\\|$', LA.norm(LA.add(x, y)), LA.norm(x) + LA.norm(y), '≤'],
          ['3. $\\|\\mathbf{x}\\| \\ge 0$', LA.norm(x), 0, '≥']
        ];

        const makeRows = (rows) => `
          <thead><tr><th>Property</th><th>LHS</th><th>RHS</th><th>Check</th></tr></thead>
          <tbody>
            ${rows.map(([p, l, r, op]) => `
              <tr>
                <td>${p}</td>
                <td><strong>${l.toFixed(2)}</strong></td>
                <td><strong>${r.toFixed(2)}</strong></td>
                <td style="color:#15803D; font-weight:700;">✓ (${op})</td>
              </tr>
            `).join('')}
          </tbody>
        `;

        const ipTable = container.querySelector('#m1IpPropTable');
        const normTable = container.querySelector('#m1NormPropTable');
        ipTable.innerHTML = makeRows(ipRows);
        normTable.innerHTML = makeRows(normRows);
        renderMath(container.querySelector('#m1PropVectorsLabel'));
        renderMath(ipTable);
        renderMath(normTable);
      }

      container.querySelector('#m1NewPropTest').addEventListener('click', runPropVerifier);

      // =====================================================================
      // SECTION E: QUICK CONCEPT CHECK INTERACTION
      // =====================================================================
      const explanations = {
        1: '⟨x, y⟩ = 2(1) + (-3)(2) + 1(-1) = 2 - 6 - 1 = -5 < 0, so the angle is obtuse.',
        2: 'By norm property 1, ‖c x‖ = |c| · ‖x‖. Since |-3| = 3, ‖-3x‖ = 3 · 5 = 15 (a norm is never negative!).',
        3: 'First find ‖v‖ = √(3² + (-4)² + 0²) = √25 = 5, then divide v by 5 to get a unit vector.'
      };

      container.querySelectorAll('.m1-quiz-item').forEach(item => {
        const qNum = item.dataset.q;
        const ans = item.dataset.ans;
        const fb = container.querySelector(`#m1Fb${qNum}`);
        item.querySelectorAll('.m1-opt-btn').forEach(btn => {
          btn.addEventListener('click', () => {
            item.querySelectorAll('.m1-opt-btn').forEach(b => b.classList.remove('correct', 'wrong'));
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
      drawCanvas();
      buildVectorInputs();
      runPropVerifier();
    }
  });
})();
