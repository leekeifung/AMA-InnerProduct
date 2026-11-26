/**
 * AMA1751 Linear Algebra — Topic 6: Inner Product
 * Module 7: General Inner Product Spaces, Orthogonal Polynomials & Fourier Series
 */
(function () {
  function clean(x, d = 3) {
    if (!Number.isFinite(x)) return 'NaN';
    if (Math.abs(x) < 1e-11) return '0';
    return String(Number(x.toFixed(d)));
  }

  /**
   * Composite Simpson integration with explicit interior split points
   * so piecewise functions (|t|, sgn(t), pulses) integrate cleanly.
   */
  function integrateSplit(fn, a, b, splits = [], panelsPerSeg = 120) {
    const pts = [a, ...splits.filter(s => s > a + 1e-12 && s < b - 1e-12).sort((x, y) => x - y), b];
    let total = 0;
    for (let seg = 0; seg < pts.length - 1; seg++) {
      const lo = pts[seg], hi = pts[seg + 1];
      const n = panelsPerSeg % 2 === 0 ? panelsPerSeg : panelsPerSeg + 1;
      const h = (hi - lo) / n;
      // Offset tiny eps inward at endpoints so jump discontinuities evaluate cleanly
      const eps = 1e-11 * (hi - lo);
      let s = fn(lo + eps) + fn(hi - eps);
      for (let i = 1; i < n; i++) {
        const x = lo + i * h;
        s += (i % 2 === 0 ? 2 : 4) * fn(x);
      }
      total += (h / 3) * s;
    }
    return total;
  }

  /** Standard Legendre polynomial P_k(t) with P_k(1) = 1 via 3-term recurrence */
  function legendreP(k, t) {
    if (k === 0) return 1;
    if (k === 1) return t;
    let p0 = 1, p1 = t;
    for (let n = 1; n < k; n++) {
      const p2 = ((2 * n + 1) * t * p1 - n * p0) / (n + 1);
      p0 = p1; p1 = p2;
    }
    return p1;
  }

  AMA1751.registerModule(7, {
    render(meta) {
      return `
        <style>
          .m7-table {
            width: 100%;
            border-collapse: collapse;
            font-size: 0.86rem;
            margin-top: 10px;
          }
          .m7-table th, .m7-table td {
            border: 1px solid var(--border);
            padding: 8px 10px;
            text-align: left;
            vertical-align: top;
          }
          .m7-table th {
            background: #F8FAFC;
            font-weight: 700;
            color: #0F172A;
          }
          .m7-canvas-wrap {
            position: relative;
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 10px;
            overflow: hidden;
          }
          .m7-canvas {
            display: block;
            width: 100%;
            height: 360px;
            touch-action: none;
          }
          .m7-canvas-sm {
            display: block;
            width: 100%;
            height: 220px;
          }
          .m7-controls {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            margin-bottom: 10px;
            align-items: center;
          }
          .m7-chip {
            padding: 6px 11px;
            border-radius: 999px;
            border: 1px solid var(--border);
            background: #F8FAFC;
            font-size: 0.82rem;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.15s;
          }
          .m7-chip:hover, .m7-chip.active {
            background: var(--polyu-red-soft);
            border-color: #E8B4C4;
            color: var(--polyu-red-dark);
          }
          .m7-panel {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 12px 14px;
            font-size: 0.88rem;
            line-height: 1.6;
          }
          .m7-slider-row {
            display: flex;
            align-items: center;
            gap: 10px;
            margin: 6px 0;
          }
          .m7-slider-row input[type="range"] { flex: 1; }
          .m7-slider-lbl {
            min-width: 112px;
            font-weight: 700;
            font-size: 0.84rem;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          }
          .m7-toggle-row {
            display: flex;
            flex-wrap: wrap;
            gap: 14px;
            align-items: center;
            margin: 6px 0 10px;
            font-size: 0.85rem;
            font-weight: 600;
            color: #334155;
          }
          .m7-toggle-row label {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            cursor: pointer;
          }
          .m7-stat-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(112px, 1fr));
            gap: 8px;
            margin: 10px 0;
          }
          .m7-stat-card {
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 8px 9px;
            text-align: center;
          }
          .m7-stat-lbl {
            font-size: 0.7rem;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.04em;
            color: var(--text-muted);
          }
          .m7-stat-val {
            font-size: 0.98rem;
            font-weight: 800;
            color: #0F172A;
            margin-top: 2px;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          }
          .m7-badge {
            display: inline-block;
            padding: 5px 10px;
            border-radius: 8px;
            font-size: 0.82rem;
            font-weight: 700;
            margin-top: 6px;
          }
          .m7-badge.ok { background: #DCFCE7; color: #166534; border: 1px solid #86EFAC; }
          .m7-badge.bad { background: #FEE2E2; color: #991B1B; border: 1px solid #FCA5A5; }
          .m7-energy-bar {
            width: 100%;
            height: 12px;
            background: #E2E8F0;
            border-radius: 999px;
            overflow: hidden;
            margin: 6px 0;
          }
          .m7-energy-fill {
            height: 100%;
            background: linear-gradient(90deg, #0F766E, #14B8A6);
            border-radius: 999px;
            transition: width 0.2s ease;
          }
          .m7-coeff-list {
            display: flex;
            flex-wrap: wrap;
            gap: 6px;
            margin-top: 8px;
          }
          .m7-coeff-pill {
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 6px;
            padding: 3px 8px;
            font-size: 0.78rem;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
            font-weight: 700;
          }
          .m7-quiz-item {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 14px 16px;
            margin-bottom: 12px;
          }
          .m7-opt-group { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 10px; }
          .m7-opt-btn {
            padding: 7px 13px;
            border: 1px solid #CBD5E1;
            border-radius: 8px;
            background: #FFFFFF;
            cursor: pointer;
            font-weight: 600;
            font-size: 0.88rem;
          }
          .m7-opt-btn:hover { background: #F1F5F9; }
          .m7-opt-btn.correct { background: #DCFCE7; border-color: #16A34A; color: #166534; }
          .m7-opt-btn.wrong { background: #FEE2E2; border-color: #DC2626; color: #991B1B; }
          .m7-fb { margin-top: 8px; font-size: 0.88rem; font-weight: 500; }
        </style>

        <!-- HEADER -->
        <div class="module-header">
          <span class="slide-ref">Part 7 · Lecture Notes ${meta.slides}</span>
          <h2>7. General Inner Product Spaces, Orthogonal Polynomials &amp; Fourier Series</h2>
          <p>Every geometric theorem from Modules 1–6—length, angle, Cauchy–Schwarz, Gram–Schmidt, and best-approximation projection—transfers <strong>unchanged</strong> to polynomials, matrices, and continuous functions once we define a valid inner product $\\langle \\mathbf{u}, \\mathbf{v} \\rangle$.</p>
        </div>

        <!-- SECTION A: AXIOMS & CLASSIC SPACES -->
        <div class="card">
          <h3>The Four Axioms &amp; Four Iconic Inner Product Spaces</h3>
          <div class="grid-2">
            <div class="callout definition">
              <strong>Definition (Real Inner Product Space):</strong><br/>
              An <strong>inner product</strong> on a real vector space $V$ assigns a scalar $\\langle \\mathbf{u}, \\mathbf{v} \\rangle \\in \\mathbb{R}$ to every pair $\\mathbf{u}, \\mathbf{v} \\in V$ satisfying:
              <ol style="margin:6px 0 0 18px; padding:0; font-size:0.9rem; line-height:1.58;">
                <li><strong>Symmetry:</strong> $\\langle \\mathbf{u}, \\mathbf{v} \\rangle = \\langle \\mathbf{v}, \\mathbf{u} \\rangle$</li>
                <li><strong>Additivity:</strong> $\\langle \\mathbf{u} + \\mathbf{v}, \\mathbf{w} \\rangle = \\langle \\mathbf{u}, \\mathbf{w} \\rangle + \\langle \\mathbf{v}, \\mathbf{w} \\rangle$</li>
                <li><strong>Homogeneity:</strong> $\\langle c\\mathbf{u}, \\mathbf{v} \\rangle = c\\langle \\mathbf{u}, \\mathbf{v} \\rangle$</li>
                <li><strong>Positivity:</strong> $\\langle \\mathbf{u}, \\mathbf{u} \\rangle \\ge 0$, and $\\langle \\mathbf{u}, \\mathbf{u} \\rangle = 0 \\iff \\mathbf{u} = \\mathbf{0}$</li>
              </ol>
            </div>

            <div class="callout tip">
              <strong>Universal Theorems (Hold in <em>Every</em> Inner Product Space!):</strong>
              <ul style="margin:6px 0 0 18px; padding:0; font-size:0.9rem; line-height:1.6;">
                <li><strong>Length &amp; Distance:</strong> $\\|\\mathbf{u}\\| = \\sqrt{\\langle \\mathbf{u}, \\mathbf{u} \\rangle}$, &nbsp;$d(\\mathbf{u}, \\mathbf{v}) = \\|\\mathbf{u} - \\mathbf{v}\\|$.</li>
                <li><strong>Cauchy–Schwarz Inequality:</strong> $|\\langle \\mathbf{u}, \\mathbf{v} \\rangle| \\le \\|\\mathbf{u}\\|\\|\\mathbf{v}\\|$.</li>
                <li><strong>Triangle Inequality:</strong> $\\|\\mathbf{u} + \\mathbf{v}\\| \\le \\|\\mathbf{u}\\| + \\|\\mathbf{v}\\|$.</li>
                <li><strong>Orthogonal Projection onto $W = \\operatorname{span}\\{\\mathbf{w}_1,\\dots,\\mathbf{w}_k\\}$</strong> (orthogonal basis):
                  $$\\operatorname{proj}_W(\\mathbf{v}) = \\sum_{i=1}^k \\frac{\\langle \\mathbf{v}, \\mathbf{w}_i \\rangle}{\\langle \\mathbf{w}_i, \\mathbf{w}_i \\rangle}\\mathbf{w}_i$$
                </li>
              </ul>
            </div>
          </div>

          <div style="overflow-x:auto;">
            <table class="m7-table">
              <thead>
                <tr>
                  <th>Vector Space $V$</th>
                  <th>Inner Product $\\langle \\mathbf{u}, \\mathbf{v} \\rangle$</th>
                  <th>Condition for Positivity</th>
                  <th>Where It Shines</th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <td><strong>Weighted $\\mathbb{R}^n$</strong></td>
                  <td>$\\mathbf{u}^T W \\mathbf{v}$ (e.g. $\\sum w_i u_i v_i$)</td>
                  <td>$W$ symmetric positive definite ($w_i > 0$ for diagonal)</td>
                  <td>Weighted least squares, Mahalanobis distance in statistics</td>
                </tr>
                <tr>
                  <td><strong>Polynomials $\\mathbb{P}_n$</strong></td>
                  <td>$\\sum_{i=0}^m p(t_i) q(t_i)$</td>
                  <td>At least $n+1$ <em>distinct</em> sample points $t_i$</td>
                  <td>Discrete polynomial fitting &amp; discrete orthogonal polynomials</td>
                </tr>
                <tr>
                  <td><strong>Functions $C[a, b]$</strong></td>
                  <td>$\\int_a^b f(t)g(t)\\,dt$</td>
                  <td>Always valid on $C[a,b]$ (if $f(t)^2 \\ge 0$ integrates to $0$, $f \\equiv 0$)</td>
                  <td>Legendre polynomials ($[-1,1]$) &amp; Fourier series ($[-\\pi,\\pi]$)</td>
                </tr>
                <tr>
                  <td><strong>Matrices $\\mathbb{R}^{m\\times n}$</strong></td>
                  <td>$\\operatorname{tr}(A^T B) = \\sum_{i,j} a_{ij} b_{ij}$</td>
                  <td>Always valid ($\\|A\\|_F^2 = \\sum a_{ij}^2$)</td>
                  <td>Low-rank matrix approximation &amp; SVD (Frobenius norm)</td>
                </tr>
              </tbody>
            </table>
          </div>
        </div>

        <!-- SECTION B: INTERACTIVE LAB 1 — WEIGHTED INNER PRODUCT GEOMETRY STUDIO -->
        <div class="card">
          <h3>Interactive Lab 1: Weighted Inner Product Geometry Studio ($\\langle \\mathbf{u}, \\mathbf{v} \\rangle_W = \\mathbf{u}^T W \\mathbf{v}$)</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            <strong>Drag</strong> the blue ($\\mathbf{u}$) and red ($\\mathbf{v}$) arrowheads. Changing the symmetric matrix $W = \\begin{bmatrix} w_1 & c \\\\ c & w_2 \\end{bmatrix}$ warps the unit "circle" $\\{\\mathbf{x} : \\|\\mathbf{x}\\|_W = 1\\}$ into an <strong>ellipse</strong> and tilts what "$90^\\circ$ orthogonal" means!
          </p>

          <div class="m7-controls">
            <strong style="font-size:0.84rem;">Weight Matrix $W$:</strong>
            <button class="m7-chip active" data-wp="euclid">Euclidean $I$</button>
            <button class="m7-chip" data-wp="diag">Weighted $\\operatorname{diag}(4, 1)$</button>
            <button class="m7-chip" data-wp="corr">Correlated $\\begin{bmatrix}2 & 1\\\\ 1 & 2\\end{bmatrix}$</button>
            <button class="m7-chip" data-wp="invalid">Invalid Indefinite (Breaks Positivity!)</button>
          </div>

          <div class="grid-2" style="align-items:start;">
            <div class="m7-canvas-wrap">
              <canvas id="m7CanvasW" class="m7-canvas" style="cursor:grab;"></canvas>
            </div>

            <div class="m7-panel">
              <div style="font-weight:700; margin-bottom:4px;">Symmetric Matrix $W = \\begin{bmatrix} w_1 & c \\\\ c & w_2 \\end{bmatrix}$:</div>
              <div class="m7-slider-row">
                <span class="m7-slider-lbl" id="m7W1Lbl">w₁ = 1.00</span>
                <input type="range" id="m7W1" min="-2" max="5" step="0.25" value="1" />
              </div>
              <div class="m7-slider-row">
                <span class="m7-slider-lbl" id="m7W2Lbl">w₂ = 1.00</span>
                <input type="range" id="m7W2" min="-2" max="5" step="0.25" value="1" />
              </div>
              <div class="m7-slider-row">
                <span class="m7-slider-lbl" id="m7WCLbl">c  = 0.00</span>
                <input type="range" id="m7WC" min="-3" max="3" step="0.25" value="0" />
              </div>

              <div id="m7WStats" style="margin-top:8px; border-top:1px solid var(--border); padding-top:8px;"></div>
            </div>
          </div>
        </div>

        <!-- SECTION C: INTERACTIVE LAB 2 — LEGENDRE ORTHOGONAL POLYNOMIALS IN C[-1, 1] -->
        <div class="card">
          <h3>Interactive Lab 2: Function Space $C[-1, 1]$ &amp; Legendre Orthogonal Polynomials</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            Applying <strong>Gram–Schmidt</strong> to $\\{1, t, t^2, t^3, \\dots\\}$ under $\\langle f, g \\rangle = \\int_{-1}^1 f(t)g(t)\\,dt$ produces the orthogonal <strong>Legendre polynomials</strong> $P_0(t)=1,\\; P_1(t)=t,\\; P_2(t)=\\frac{1}{2}(3t^2-1),\\; P_3(t)=\\frac{1}{2}(5t^3-3t),\\dots$ with $\\|P_k\\|^2 = \\frac{2}{2k+1}$. Project any function onto $\\operatorname{span}\\{P_0,\\dots,P_N\\}$ below!
          </p>

          <div class="m7-controls">
            <strong style="font-size:0.84rem;">Target $f(t)$:</strong>
            <button class="m7-chip active" data-lf="exp">$e^t$</button>
            <button class="m7-chip" data-lf="abs">$|t|$ (Corner)</button>
            <button class="m7-chip" data-lf="sin">$\\sin(\\pi t)$</button>
            <button class="m7-chip" data-lf="sgn">$\\operatorname{sgn}(t)$ (Step)</button>
            <button class="m7-chip" data-lf="runge">$\\frac{1}{1+10t^2}$ (Bell)</button>
          </div>

          <div class="m7-slider-row" style="max-width:520px;">
            <span class="m7-slider-lbl" id="m7DegLbl">Degree N = 2</span>
            <input type="range" id="m7Deg" min="0" max="6" step="1" value="2" />
          </div>

          <div class="m7-toggle-row">
            <label><input type="checkbox" id="m7ShowTaylor" checked /> Compare with Maclaurin/Taylor Polynomial (when analytic)</label>
            <label><input type="checkbox" id="m7ShowBasis" /> Show Basis Curves $P_0(t) \\dots P_4(t)$ Only</label>
          </div>

          <div class="grid-2" style="align-items:start;">
            <div class="m7-canvas-wrap">
              <canvas id="m7CanvasLeg" class="m7-canvas"></canvas>
            </div>

            <div class="m7-panel" id="m7LegStats"></div>
          </div>
        </div>

        <!-- SECTION D: INTERACTIVE LAB 3 — FOURIER SERIES AS ORTHOGONAL PROJECTION -->
        <div class="card">
          <h3>Interactive Lab 3: Fourier Approximation Studio in $C[-\\pi, \\pi]$</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            Under $\\langle f, g \\rangle = \\int_{-\\pi}^{\\pi} f(t)g(t)\\,dt$, the trigonometric functions $\\{1, \\cos t, \\sin t, \\cos 2t, \\sin 2t, \\dots\\}$ are <strong>mutually orthogonal</strong>! Therefore, the $N$-th order <strong>Fourier series</strong> $S_N(t)$ is literally the orthogonal projection $\\operatorname{proj}_{W_N}(f)$:
            $$S_N(t) = \\frac{a_0}{2} + \\sum_{k=1}^N \\bigl(a_k \\cos(kt) + b_k \\sin(kt)\\bigr), \\qquad a_k = \\frac{1}{\\pi}\\int_{-\\pi}^{\\pi} f(t)\\cos(kt)\\,dt, \\quad b_k = \\frac{1}{\\pi}\\int_{-\\pi}^{\\pi} f(t)\\sin(kt)\\,dt$$
          </p>

          <div class="m7-controls">
            <strong style="font-size:0.84rem;">Periodic Signal $f(t)$:</strong>
            <button class="m7-chip active" data-ff="square">Square Wave (Gibbs!)</button>
            <button class="m7-chip" data-ff="sawtooth">Sawtooth $t/\\pi$</button>
            <button class="m7-chip" data-ff="triangle">Triangle Wave</button>
            <button class="m7-chip" data-ff="rectified">Half-Wave Rectified $\\sin t$</button>
            <button class="m7-chip" data-ff="pulse">Symmetric Pulse</button>
            <button class="m7-chip" id="m7AnimBtn" style="margin-left:auto; background:#CCFBF1; border-color:#5EEAD4; color:#115E59;">▶ Animate $N$</button>
          </div>

          <div class="m7-slider-row" style="max-width:560px;">
            <span class="m7-slider-lbl" id="m7HarmLbl">Harmonics N = 3</span>
            <input type="range" id="m7Harm" min="1" max="25" step="1" value="3" />
          </div>

          <div class="grid-2" style="align-items:start;">
            <div class="m7-canvas-wrap">
              <canvas id="m7CanvasFour" class="m7-canvas"></canvas>
            </div>

            <div class="m7-panel">
              <div style="font-weight:700; margin-bottom:4px;">Harmonic Spectrum ($|a_k|$ teal, $|b_k|$ crimson):</div>
              <div class="m7-canvas-wrap" style="margin-bottom:10px;">
                <canvas id="m7CanvasSpec" class="m7-canvas-sm"></canvas>
              </div>
              <div id="m7FourStats"></div>
            </div>
          </div>

          <details class="code-details">
            <summary>Python (SymPy &amp; NumPy) Companion Code — Symbolic Gram–Schmidt for Legendre Polynomials &amp; Fourier Projection</summary>
            <pre><code>import sympy as sp

t = sp.Symbol('t', real=True)

# 1) Exact Inner Product on C[-1, 1] and Gram-Schmidt -> Legendre Polynomials
def ip_L2(f, g):
    return sp.integrate(f * g, (t, -1, 1))

monomials = [1, t, t**2, t**3]
orth = []
for v in monomials:
    w = v
    for u in orth:
        w = w - (ip_L2(v, u) / ip_L2(u, u)) * u
    orth.append(sp.simplify(w))

# Standard Legendre scaling P_k(1) = 1
P = [sp.simplify(w / w.subs(t, 1)) for w in orth]
for k, pk in enumerate(P):
    print(f"P_{k}(t) = {pk},  ||P_{k}||^2 = {ip_L2(pk, pk)}")

# 2) Best L2 Degree-3 Projection of f(t) = exp(t) on [-1, 1]
f = sp.exp(t)
f_hat = sum((ip_L2(f, pk) / ip_L2(pk, pk)) * pk for pk in P)
print("L2 projection of e^t (deg 3):", sp.N(sp.expand(f_hat), 5))

# 3) Exact Fourier Projection of Square Wave sgn(t) on [-pi, pi]
k = sp.Symbol('k', integer=True, positive=True)
b_k = (1 / sp.pi) * (
    sp.integrate((-1) * sp.sin(k * t), (t, -sp.pi, 0)) +
    sp.integrate((+1) * sp.sin(k * t), (t, 0, sp.pi))
)
print("Square wave b_k =", sp.simplify(b_k))  # 2*(1 - (-1)**k) / (pi*k)</code></pre>
          </details>
        </div>

        <!-- SECTION E: QUICK CONCEPT CHECK -->
        <div class="card">
          <h3>Quick Concept Check (Exam Pitfall Radar)</h3>

          <div class="m7-quiz-item" data-q="1" data-ans="C">
            <strong>Q1.</strong> On the polynomial space $\\mathbb{P}_2$ (polynomials of degree $\\le 2$), consider the formula $\\langle p, q \\rangle = p(0)q(0) + p(1)q(1)$. Does this define a valid inner product on $\\mathbb{P}_2$?
            <div class="m7-opt-group">
              <button class="m7-opt-btn" data-opt="A">Yes — it satisfies all four axioms.</button>
              <button class="m7-opt-btn" data-opt="B">No — it fails the Symmetry axiom.</button>
              <button class="m7-opt-btn" data-opt="C">No — it fails Positivity: the nonzero polynomial $p(t) = t(t-1) \\in \\mathbb{P}_2$ has $\\langle p, p \\rangle = 0$.</button>
              <button class="m7-opt-btn" data-opt="D">No — it fails Additivity because polynomials are nonlinear.</button>
            </div>
            <div class="m7-fb" id="m7Fb1"></div>
          </div>

          <div class="m7-quiz-item" data-q="2" data-ans="B">
            <strong>Q2.</strong> On $\\mathbb{R}^2$ equipped with the weighted inner product $\\langle \\mathbf{u}, \\mathbf{v} \\rangle = 4u_1 v_1 + u_2 v_2$, which vector is <em>orthogonal</em> to $\\mathbf{u} = (1, 2)^T$?
            <div class="m7-opt-group">
              <button class="m7-opt-btn" data-opt="A">$(-2, 1)^T$</button>
              <button class="m7-opt-btn" data-opt="B">$(1, -2)^T$</button>
              <button class="m7-opt-btn" data-opt="C">$(2, -1)^T$</button>
              <button class="m7-opt-btn" data-opt="D">$(1, -4)^T$</button>
            </div>
            <div class="m7-fb" id="m7Fb2"></div>
          </div>

          <div class="m7-quiz-item" data-q="3" data-ans="A">
            <strong>Q3.</strong> You already found the best degree-2 $L^2$ approximation $\\hat{f}_2 = c_0 P_0 + c_1 P_1 + c_2 P_2$ on $[-1, 1]$ using the orthogonal Legendre polynomials. Now you upgrade to the best degree-3 approximation $\\hat{f}_3$. What happens to $c_0, c_1, c_2$?
            <div class="m7-opt-group">
              <button class="m7-opt-btn" data-opt="A">They remain completely unchanged—just add $c_3 P_3$, where $c_3 = \\frac{\\langle f, P_3 \\rangle}{\\langle P_3, P_3 \\rangle}$.</button>
              <button class="m7-opt-btn" data-opt="B">All four coefficients must be recomputed by solving a new $4\\times 4$ linear system.</button>
              <button class="m7-opt-btn" data-opt="C">Only $c_0$ stays unchanged; $c_1$ and $c_2$ must be rescaled.</button>
              <button class="m7-opt-btn" data-opt="D">The coefficients $c_0, c_1, c_2$ are halved.</button>
            </div>
            <div class="m7-fb" id="m7Fb3"></div>
          </div>

          <div class="m7-quiz-item" data-q="4" data-ans="D">
            <strong>Q4.</strong> Applying the Cauchy–Schwarz inequality $|\\langle f, g \\rangle| \\le \\|f\\|\\|g\\|$ to a continuous function $f$ and the constant function $g(t) = 1$ in $C[0, 1]$ with $\\langle f, g \\rangle = \\int_0^1 f(t)g(t)\\,dt$ guarantees which inequality?
            <div class="m7-opt-group">
              <button class="m7-opt-btn" data-opt="A">$\\int_0^1 f(t)^2\\,dt \\le \\left(\\int_0^1 f(t)\\,dt\\right)^2$</button>
              <button class="m7-opt-btn" data-opt="B">$\\int_0^1 f(t)\\,dt \\le 1$</button>
              <button class="m7-opt-btn" data-opt="C">$\\left(\\int_0^1 f(t)\\,dt\\right)^2 = \\int_0^1 f(t)^2\\,dt$</button>
              <button class="m7-opt-btn" data-opt="D">$\\left(\\int_0^1 f(t)\\,dt\\right)^2 \\le \\int_0^1 f(t)^2\\,dt$</button>
            </div>
            <div class="m7-fb" id="m7Fb4"></div>
          </div>
        </div>
      `;
    },

    mount(container, { renderMath }) {
      // =====================================================================
      // LAB 1: WEIGHTED INNER PRODUCT GEOMETRY STUDIO
      // =====================================================================
      const cW = container.querySelector('#m7CanvasW');
      const ctxW = cW.getContext('2d');
      let w1 = 1, w2 = 1, wc = 0;
      let vecU = [1.4, 0.7], vecV = [-0.6, 1.5];
      let dragVec = null;

      const wPresets = {
        euclid:  { w1: 1, w2: 1, wc: 0,    u: [1.4, 0.7],  v: [-0.7, 1.4] },
        diag:    { w1: 4, w2: 1, wc: 0,    u: [0.8, 1.2],  v: [0.7, -1.5] },
        corr:    { w1: 2, w2: 2, wc: 1,    u: [1.2, 0.4],  v: [-0.9, 1.2] },
        invalid: { w1: 1, w2: -1, wc: 0.5, u: [1.2, 0.8],  v: [-0.8, 1.4] }
      };

      function ipW(a, b) {
        return w1 * a[0] * b[0] + wc * (a[0] * b[1] + a[1] * b[0]) + w2 * a[1] * b[1];
      }

      function layoutW() {
        const r = cW.getBoundingClientRect();
        const w = r.width || 440, h = r.height || 360;
        const sc = Math.min(w, h) / 6.2;
        return { w, h, ox: w / 2, oy: h / 2, sc };
      }
      function toSW(x, y, L) { return [L.ox + x * L.sc, L.oy - y * L.sc]; }
      function toWW(px, py, L) {
        return [
          Math.max(-2.8, Math.min(2.8, (px - L.ox) / L.sc)),
          Math.max(-2.8, Math.min(2.8, (L.oy - py) / L.sc))
        ];
      }

      function drawArrowW(ctx, x1, y1, x2, y2, color, width = 3, dash = []) {
        ctx.save();
        ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = width;
        ctx.setLineDash(dash);
        ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); ctx.stroke();
        ctx.setLineDash([]);
        const ang = Math.atan2(y2 - y1, x2 - x1), hl = 10;
        if (Math.hypot(x2 - x1, y2 - y1) > 6) {
          ctx.beginPath();
          ctx.moveTo(x2, y2);
          ctx.lineTo(x2 - hl * Math.cos(ang - 0.36), y2 - hl * Math.sin(ang - 0.36));
          ctx.lineTo(x2 - hl * Math.cos(ang + 0.36), y2 - hl * Math.sin(ang + 0.36));
          ctx.closePath(); ctx.fill();
        }
        ctx.restore();
      }

      function syncWSliders() {
        container.querySelector('#m7W1').value = w1;
        container.querySelector('#m7W2').value = w2;
        container.querySelector('#m7WC').value = wc;
        container.querySelector('#m7W1Lbl').textContent = `w₁ = ${w1.toFixed(2)}`;
        container.querySelector('#m7W2Lbl').textContent = `w₂ = ${w2.toFixed(2)}`;
        container.querySelector('#m7WCLbl').textContent = `c  = ${wc.toFixed(2)}`;
      }

      function drawWLab() {
        if (!cW.isConnected) return;
        const dpr = window.devicePixelRatio || 1;
        const L = layoutW();
        cW.width = Math.round(L.w * dpr);
        cW.height = Math.round(L.h * dpr);
        ctxW.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctxW.clearRect(0, 0, L.w, L.h);

        // Grid & axes
        ctxW.strokeStyle = '#F1F5F9'; ctxW.lineWidth = 1;
        for (let g = -3; g <= 3; g++) {
          const [gx] = toSW(g, 0, L), [, gy] = toSW(0, g, L);
          ctxW.beginPath(); ctxW.moveTo(gx, 0); ctxW.lineTo(gx, L.h); ctxW.stroke();
          ctxW.beginPath(); ctxW.moveTo(0, gy); ctxW.lineTo(L.w, gy); ctxW.stroke();
        }
        ctxW.strokeStyle = '#CBD5E1'; ctxW.lineWidth = 1.4;
        ctxW.beginPath(); ctxW.moveTo(0, L.oy); ctxW.lineTo(L.w, L.oy); ctxW.stroke();
        ctxW.beginPath(); ctxW.moveTo(L.ox, 0); ctxW.lineTo(L.ox, L.h); ctxW.stroke();

        // Reference Euclidean unit circle (faint gray dashed)
        ctxW.save();
        ctxW.strokeStyle = '#CBD5E1'; ctxW.lineWidth = 1.3; ctxW.setLineDash([4, 4]);
        ctxW.beginPath(); ctxW.arc(L.ox, L.oy, L.sc, 0, Math.PI * 2); ctxW.stroke();
        ctxW.restore();

        const detW = w1 * w2 - wc * wc;
        const isSPD = w1 > 1e-6 && detW > 1e-6;
        const statsEl = container.querySelector('#m7WStats');

        if (isSPD) {
          // Draw weighted unit "circle" (ellipse): x^T W x = 1
          ctxW.save();
          ctxW.fillStyle = 'rgba(15, 118, 110, 0.08)';
          ctxW.strokeStyle = '#0F766E'; ctxW.lineWidth = 2.2;
          ctxW.beginPath();
          const steps = 160;
          for (let i = 0; i <= steps; i++) {
            const th = (i / steps) * Math.PI * 2;
            const dx = Math.cos(th), dy = Math.sin(th);
            const q = w1 * dx * dx + 2 * wc * dx * dy + w2 * dy * dy;
            const r = 1 / Math.sqrt(q);
            const [sx, sy] = toSW(r * dx, r * dy, L);
            if (i === 0) ctxW.moveTo(sx, sy); else ctxW.lineTo(sx, sy);
          }
          ctxW.closePath(); ctxW.fill(); ctxW.stroke();
          ctxW.restore();

          // Draw W-orthogonal complement line u^{perp_W}: (W u)^T x = 0
          const wu0 = w1 * vecU[0] + wc * vecU[1];
          const wu1 = wc * vecU[0] + w2 * vecU[1];
          const wuLen = Math.hypot(wu0, wu1);
          if (wuLen > 1e-6) {
            const px = -wu1 / wuLen, py = wu0 / wuLen;
            const [lx1, ly1] = toSW(-4 * px, -4 * py, L);
            const [lx2, ly2] = toSW(4 * px, 4 * py, L);
            ctxW.save();
            ctxW.strokeStyle = 'rgba(217, 119, 6, 0.72)'; ctxW.lineWidth = 1.8; ctxW.setLineDash([6, 4]);
            ctxW.beginPath(); ctxW.moveTo(lx1, ly1); ctxW.lineTo(lx2, ly2); ctxW.stroke();
            ctxW.fillStyle = '#B45309'; ctxW.font = 'bold 11px sans-serif';
            const [tx, ty] = toSW(2.15 * px, 2.15 * py, L);
            ctxW.fillText('u⊥(W)', tx + 4, ty - 4);
            ctxW.restore();
          }

          // W-orthogonal projection of v onto u
          const uu = ipW(vecU, vecU);
          const uv = ipW(vecU, vecV);
          const vv = ipW(vecV, vecV);
          const proj = uu > 1e-9 ? [vecU[0] * (uv / uu), vecU[1] * (uv / uu)] : [0, 0];

          const [pSx, pSy] = toSW(proj[0], proj[1], L);
          const [vSx, vSy] = toSW(vecV[0], vecV[1], L);
          ctxW.save();
          ctxW.strokeStyle = '#D97706'; ctxW.lineWidth = 2; ctxW.setLineDash([4, 3]);
          ctxW.beginPath(); ctxW.moveTo(pSx, pSy); ctxW.lineTo(vSx, vSy); ctxW.stroke();
          ctxW.restore();

          drawArrowW(ctxW, L.ox, L.oy, pSx, pSy, '#0F766E', 3.2);

          const nu = Math.sqrt(uu), nv = Math.sqrt(vv);
          const cosW = (nu > 1e-9 && nv > 1e-9) ? Math.max(-1, Math.min(1, uv / (nu * nv))) : 0;
          const degW = (Math.acos(cosW) * 180) / Math.PI;
          const euclDot = vecU[0] * vecV[0] + vecU[1] * vecV[1];
          const euclCos = (Math.hypot(...vecU) * Math.hypot(...vecV)) > 1e-9
            ? Math.max(-1, Math.min(1, euclDot / (Math.hypot(...vecU) * Math.hypot(...vecV)))) : 0;
          const degE = (Math.acos(euclCos) * 180) / Math.PI;

          statsEl.innerHTML = `
            <div class="m7-badge ok">✓ Valid Inner Product ($w_1 > 0$ and $\\det W = ${clean(detW, 2)} > 0$)</div>
            <div class="m7-stat-grid">
              <div class="m7-stat-card"><div class="m7-stat-lbl">$\\langle \\mathbf{u}, \\mathbf{v} \\rangle_W$</div><div class="m7-stat-val">${clean(uv, 2)}</div></div>
              <div class="m7-stat-card"><div class="m7-stat-lbl">$\\|\\mathbf{u}\\|_W,\\; \\|\\mathbf{v}\\|_W$</div><div class="m7-stat-val">${clean(nu, 2)}, ${clean(nv, 2)}</div></div>
              <div class="m7-stat-card"><div class="m7-stat-lbl">$W$-Angle $\\theta_W$</div><div class="m7-stat-val" style="color:${Math.abs(degW - 90) < 3 ? '#15803D' : '#0F172A'};">${clean(degW, 1)}°</div></div>
              <div class="m7-stat-card"><div class="m7-stat-lbl">Euclidean Angle</div><div class="m7-stat-val">${clean(degE, 1)}°</div></div>
            </div>
            <div style="font-size:0.82rem; color:#334155;">
              <strong>Cauchy–Schwarz:</strong> $|\\langle \\mathbf{u},\\mathbf{v}\\rangle_W| = ${clean(Math.abs(uv), 2)} \\le \\|\\mathbf{u}\\|_W\\|\\mathbf{v}\\|_W = ${clean(nu * nv, 2)}$ ✓<br/>
              <strong>$W$-Projection (teal):</strong> $\\operatorname{proj}_{\\mathbf{u}}^W(\\mathbf{v}) = (${clean(proj[0], 2)},\\; ${clean(proj[1], 2)})^T$ — notice the amber residual is parallel to $\\mathbf{u}^{\\perp_W}$, not Euclidean $90^\\circ$!
            </div>
          `;
        } else {
          // Construct explicit witness vector z != 0 with z^T W z <= 0
          let z;
          if (w1 <= 0) z = [1, 0];
          else if (w2 <= 0) z = [0, 1];
          else z = [-wc, w1];
          const zNorm = Math.hypot(z[0], z[1]);
          const zPlot = [1.6 * z[0] / zNorm, 1.6 * z[1] / zNorm];
          const zVal = ipW(z, z);

          const [zx, zy] = toSW(zPlot[0], zPlot[1], L);
          drawArrowW(ctxW, L.ox, L.oy, zx, zy, '#991B1B', 3, [5, 3]);
          ctxW.fillStyle = '#991B1B'; ctxW.font = 'bold 12px sans-serif';
          ctxW.fillText('Witness z (⟨z,z⟩ ≤ 0)', zx + 6, zy - 6);

          statsEl.innerHTML = `
            <div class="m7-badge bad">✗ NOT an Inner Product! Positivity Axiom Violated ($\\det W = ${clean(detW, 2)}$)</div>
            <div style="margin-top:8px; font-size:0.84rem; color:#7F1D1D;">
              <strong>Counterexample Witness:</strong> Take the nonzero vector $\\mathbf{z} = (${clean(z[0], 2)},\\; ${clean(z[1], 2)})^T \\ne \\mathbf{0}$.<br/>
              Then $\\langle \\mathbf{z}, \\mathbf{z} \\rangle_W = \\mathbf{z}^T W \\mathbf{z} = ${clean(zVal, 3)} \\le 0$!
              Length $\\|\\mathbf{z}\\|_W = \\sqrt{\\langle \\mathbf{z}, \\mathbf{z}\\rangle_W}$ is zero or imaginary, so $W$ fails Axiom 4.
            </div>
          `;
        }

        // Draw draggable u and v on top
        const [ux, uy] = toSW(vecU[0], vecU[1], L);
        const [vx, vy] = toSW(vecV[0], vecV[1], L);
        drawArrowW(ctxW, L.ox, L.oy, ux, uy, '#1D4ED8', 3.2);
        drawArrowW(ctxW, L.ox, L.oy, vx, vy, '#DC2626', 3.2);

        ctxW.fillStyle = '#1D4ED8'; ctxW.font = 'bold 13px sans-serif';
        ctxW.fillText('u', ux + 8, uy - 6);
        ctxW.fillStyle = '#DC2626';
        ctxW.fillText('v', vx + 8, vy - 6);

        ctxW.fillStyle = 'rgba(15, 23, 42, 0.55)'; ctxW.font = '11px sans-serif';
        ctxW.fillText('Drag u (blue) or v (red) · Teal ellipse = {x : ‖x‖_W = 1}', 10, 16);

        renderMath(statsEl);
      }

      function ptrW(e) {
        const r = cW.getBoundingClientRect();
        const cx = e.touches ? e.touches[0].clientX : e.clientX;
        const cy = e.touches ? e.touches[0].clientY : e.clientY;
        return [cx - r.left, cy - r.top];
      }
      function startDragW(e) {
        const L = layoutW(), [px, py] = ptrW(e);
        const [ux, uy] = toSW(vecU[0], vecU[1], L);
        const [vx, vy] = toSW(vecV[0], vecV[1], L);
        const du = Math.hypot(px - ux, py - uy);
        const dv = Math.hypot(px - vx, py - vy);
        if (Math.min(du, dv) < 22) {
          dragVec = du <= dv ? 'u' : 'v';
          if (e.touches) e.preventDefault();
        }
      }
      function moveDragW(e) {
        if (!dragVec || !cW.isConnected) return;
        if (e.touches) e.preventDefault();
        const L = layoutW(), [px, py] = ptrW(e);
        const [wx, wy] = toWW(px, py, L);
        if (dragVec === 'u') vecU = [Math.round(wx * 20) / 20, Math.round(wy * 20) / 20];
        else vecV = [Math.round(wx * 20) / 20, Math.round(wy * 20) / 20];
        drawWLab();
      }
      cW.addEventListener('mousedown', startDragW);
      cW.addEventListener('touchstart', startDragW, { passive: false });
      window.addEventListener('mousemove', moveDragW);
      window.addEventListener('touchmove', moveDragW, { passive: false });
      window.addEventListener('mouseup', () => { dragVec = null; });
      window.addEventListener('touchend', () => { dragVec = null; });

      container.querySelector('#m7W1').addEventListener('input', e => { w1 = parseFloat(e.target.value); syncWSliders(); drawWLab(); });
      container.querySelector('#m7W2').addEventListener('input', e => { w2 = parseFloat(e.target.value); syncWSliders(); drawWLab(); });
      container.querySelector('#m7WC').addEventListener('input', e => { wc = parseFloat(e.target.value); syncWSliders(); drawWLab(); });
      container.querySelectorAll('[data-wp]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-wp]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          const p = wPresets[btn.dataset.wp];
          w1 = p.w1; w2 = p.w2; wc = p.wc;
          vecU = p.u.slice(); vecV = p.v.slice();
          syncWSliders(); drawWLab();
        });
      });

      // =====================================================================
      // LAB 2: LEGENDRE ORTHOGONAL POLYNOMIALS IN C[-1, 1]
      // =====================================================================
      const cLeg = container.querySelector('#m7CanvasLeg');
      const ctxLeg = cLeg.getContext('2d');

      const legFuncs = {
        exp: {
          label: 'e^t',
          fn: t => Math.exp(t),
          splits: [],
          yRange: [-0.4, 3.0],
          taylor: (t, deg) => {
            let s = 1, term = 1;
            for (let k = 1; k <= deg; k++) { term *= t / k; s += term; }
            return s;
          }
        },
        abs: {
          label: '|t|',
          fn: t => Math.abs(t),
          splits: [0],
          yRange: [-0.25, 1.35],
          taylor: null
        },
        sin: {
          label: '\\sin(\\pi t)',
          fn: t => Math.sin(Math.PI * t),
          splits: [],
          yRange: [-1.45, 1.45],
          taylor: (t, deg) => {
            let s = 0, fact = 1;
            for (let k = 1; k <= deg; k++) {
              fact *= k;
              if (k % 2 === 1) {
                const sign = ((k - 1) / 2) % 2 === 0 ? 1 : -1;
                s += sign * Math.pow(Math.PI * t, k) / fact;
              }
            }
            return s;
          }
        },
        sgn: {
          label: '\\operatorname{sgn}(t)',
          fn: t => (t >= 0 ? 1 : -1),
          splits: [0],
          yRange: [-1.45, 1.45],
          taylor: null
        },
        runge: {
          label: '\\frac{1}{1+10t^2}',
          fn: t => 1 / (1 + 10 * t * t),
          splits: [],
          yRange: [-0.35, 1.35],
          taylor: null
        }
      };

      let activeLf = 'exp';
      let legDeg = 2;
      let showTaylor = true;
      let showBasisOnly = false;

      function drawLegLab() {
        if (!cLeg.isConnected) return;
        const dpr = window.devicePixelRatio || 1;
        const r = cLeg.getBoundingClientRect();
        const w = r.width || 440, h = r.height || 360;
        cLeg.width = Math.round(w * dpr);
        cLeg.height = Math.round(h * dpr);
        ctxLeg.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctxLeg.clearRect(0, 0, w, h);

        const padL = 40, padR = 16, padT = 24, padB = 30;
        const plotW = w - padL - padR, plotH = h - padT - padB;
        const lf = legFuncs[activeLf];
        const [yMin, yMax] = showBasisOnly ? [-1.25, 1.25] : lf.yRange;

        const mapX = t => padL + ((t + 1) / 2) * plotW;
        const mapY = y => padT + ((yMax - y) / (yMax - yMin)) * plotH;

        // Grid & axes
        ctxLeg.strokeStyle = '#F1F5F9'; ctxLeg.lineWidth = 1;
        [-1, -0.5, 0, 0.5, 1].forEach(t => {
          const x = mapX(t);
          ctxLeg.beginPath(); ctxLeg.moveTo(x, padT); ctxLeg.lineTo(x, h - padB); ctxLeg.stroke();
        });
        ctxLeg.strokeStyle = '#94A3B8'; ctxLeg.lineWidth = 1.3;
        if (yMin <= 0 && yMax >= 0) {
          ctxLeg.beginPath(); ctxLeg.moveTo(padL, mapY(0)); ctxLeg.lineTo(w - padR, mapY(0)); ctxLeg.stroke();
        }
        ctxLeg.beginPath(); ctxLeg.moveTo(mapX(0), padT); ctxLeg.lineTo(mapX(0), h - padB); ctxLeg.stroke();

        // Axis labels
        ctxLeg.fillStyle = '#64748B'; ctxLeg.font = '11px sans-serif';
        [-1, -0.5, 0, 0.5, 1].forEach(t => ctxLeg.fillText(String(t), mapX(t) - 7, h - 10));

        const statsEl = container.querySelector('#m7LegStats');

        if (showBasisOnly) {
          const colors = ['#1D4ED8', '#0F766E', '#D97706', '#DC2626', '#7E22CE'];
          for (let k = 0; k <= 4; k++) {
            ctxLeg.save();
            ctxLeg.strokeStyle = colors[k]; ctxLeg.lineWidth = 2.4;
            ctxLeg.beginPath();
            for (let i = 0; i <= 180; i++) {
              const t = -1 + (2 * i) / 180;
              const x = mapX(t), y = mapY(legendreP(k, t));
              if (i === 0) ctxLeg.moveTo(x, y); else ctxLeg.lineTo(x, y);
            }
            ctxLeg.stroke();
            ctxLeg.restore();
          }
          ctxLeg.fillStyle = '#0F172A'; ctxLeg.font = 'bold 11px sans-serif';
          ctxLeg.fillText('Legendre Orthogonal Basis P₀ (blue), P₁ (teal), P₂ (amber), P₃ (red), P₄ (purple)', padL, 15);

          statsEl.innerHTML = `
            <div><strong>First Five Legendre Polynomials (Gram–Schmidt on $\\{1, t, t^2, t^3, t^4\\}$):</strong></div>
            <ul style="margin:6px 0 8px 18px; padding:0; font-size:0.86rem; line-height:1.65;">
              <li><span style="color:#1D4ED8; font-weight:700;">$P_0(t) = 1$</span>, &nbsp;$\\|P_0\\|^2 = 2$</li>
              <li><span style="color:#0F766E; font-weight:700;">$P_1(t) = t$</span>, &nbsp;$\\|P_1\\|^2 = \\frac{2}{3}$</li>
              <li><span style="color:#D97706; font-weight:700;">$P_2(t) = \\frac{1}{2}(3t^2 - 1)$</span>, &nbsp;$\\|P_2\\|^2 = \\frac{2}{5}$</li>
              <li><span style="color:#DC2626; font-weight:700;">$P_3(t) = \\frac{1}{2}(5t^3 - 3t)$</span>, &nbsp;$\\|P_3\\|^2 = \\frac{2}{7}$</li>
              <li><span style="color:#7E22CE; font-weight:700;">$P_4(t) = \\frac{1}{8}(35t^4 - 30t^2 + 3)$</span>, &nbsp;$\\|P_4\\|^2 = \\frac{2}{9}$</li>
            </ul>
            <div class="m7-badge ok">✓ Mutual Orthogonality: $\\int_{-1}^1 P_i(t) P_j(t)\\,dt = 0$ for all $i \\ne j$</div>
          `;
          renderMath(statsEl);
          return;
        }

        // Compute Legendre projection coefficients c_k = ((2k+1)/2) * int_{-1}^1 f(t) P_k(t) dt
        const coeffs = [];
        for (let k = 0; k <= legDeg; k++) {
          const ip = integrateSplit(t => lf.fn(t) * legendreP(k, t), -1, 1, lf.splits, 100);
          coeffs.push(((2 * k + 1) / 2) * ip);
        }
        const evalProj = t => {
          let s = 0;
          for (let k = 0; k <= legDeg; k++) s += coeffs[k] * legendreP(k, t);
          return s;
        };

        // Shaded L2 error between f(t) and proj(t)
        ctxLeg.save();
        ctxLeg.fillStyle = 'rgba(217, 119, 6, 0.14)';
        ctxLeg.beginPath();
        const samples = 220;
        for (let i = 0; i <= samples; i++) {
          const t = -1 + (2 * i) / samples;
          const x = mapX(t), y = mapY(lf.fn(t));
          if (i === 0) ctxLeg.moveTo(x, y); else ctxLeg.lineTo(x, y);
        }
        for (let i = samples; i >= 0; i--) {
          const t = -1 + (2 * i) / samples;
          ctxLeg.lineTo(mapX(t), mapY(evalProj(t)));
        }
        ctxLeg.closePath(); ctxLeg.fill();
        ctxLeg.restore();

        // Draw Taylor polynomial if enabled & available
        let taylorErrSq = null;
        if (showTaylor && lf.taylor) {
          taylorErrSq = integrateSplit(t => Math.pow(lf.fn(t) - lf.taylor(t, legDeg), 2), -1, 1, lf.splits, 100);
          ctxLeg.save();
          ctxLeg.strokeStyle = '#7E22CE'; ctxLeg.lineWidth = 2.2; ctxLeg.setLineDash([6, 4]);
          ctxLeg.beginPath();
          for (let i = 0; i <= samples; i++) {
            const t = -1 + (2 * i) / samples;
            const yVal = Math.max(yMin - 0.5, Math.min(yMax + 0.5, lf.taylor(t, legDeg)));
            const x = mapX(t), y = mapY(yVal);
            if (i === 0) ctxLeg.moveTo(x, y); else ctxLeg.lineTo(x, y);
          }
          ctxLeg.stroke();
          ctxLeg.restore();
        }

        // Draw target f(t) (blue)
        ctxLeg.save();
        ctxLeg.strokeStyle = '#1D4ED8'; ctxLeg.lineWidth = 2.6;
        ctxLeg.beginPath();
        for (let i = 0; i <= samples; i++) {
          const t = -1 + (2 * i) / samples;
          const x = mapX(t), y = mapY(lf.fn(t));
          if (i === 0 || (activeLf === 'sgn' && Math.abs(t) < 0.012)) ctxLeg.moveTo(x, y);
          else ctxLeg.lineTo(x, y);
        }
        ctxLeg.stroke();
        ctxLeg.restore();

        // Draw L2 orthogonal projection (teal)
        ctxLeg.save();
        ctxLeg.strokeStyle = '#0F766E'; ctxLeg.lineWidth = 3;
        ctxLeg.beginPath();
        for (let i = 0; i <= samples; i++) {
          const t = -1 + (2 * i) / samples;
          const x = mapX(t), y = mapY(evalProj(t));
          if (i === 0) ctxLeg.moveTo(x, y); else ctxLeg.lineTo(x, y);
        }
        ctxLeg.stroke();
        ctxLeg.restore();

        // Legend inside canvas
        ctxLeg.fillStyle = '#0F172A'; ctxLeg.font = 'bold 11px sans-serif';
        ctxLeg.fillText(`Blue: f(t) · Teal: Best L² Projection (deg ${legDeg})${showTaylor && lf.taylor ? ' · Purple dashed: Taylor' : ''}`, padL, 15);

        // Compute exact L2 norms & Pythagorean decomposition
        const normFSq = integrateSplit(t => Math.pow(lf.fn(t), 2), -1, 1, lf.splits, 100);
        let projNormSq = 0;
        for (let k = 0; k <= legDeg; k++) projNormSq += coeffs[k] * coeffs[k] * (2 / (2 * k + 1));
        const errSq = Math.max(0, normFSq - projNormSq);
        const pct = normFSq > 1e-12 ? Math.min(100, (projNormSq / normFSq) * 100) : 100;

        const pills = coeffs.map((c, k) => `<span class="m7-coeff-pill">c<sub>${k}</sub> = ${clean(c, 3)}</span>`).join('');
        statsEl.innerHTML = `
          <div><strong>Orthogonal Projection of $f(t) = ${lf.label}$ onto $\\mathbb{P}_${legDeg}$:</strong></div>
          <div class="m7-coeff-list">${pills}</div>
          <div class="m7-stat-grid">
            <div class="m7-stat-card"><div class="m7-stat-lbl">$\\|f\\|_{L^2}^2$</div><div class="m7-stat-val">${clean(normFSq, 4)}</div></div>
            <div class="m7-stat-card"><div class="m7-stat-lbl">$\\|\\hat{f}_N\\|_{L^2}^2$</div><div class="m7-stat-val" style="color:#0F766E;">${clean(projNormSq, 4)}</div></div>
            <div class="m7-stat-card"><div class="m7-stat-lbl">Error $\\|f - \\hat{f}_N\\|^2$</div><div class="m7-stat-val" style="color:#D97706;">${clean(errSq, 4)}</div></div>
          </div>
          <div style="font-size:0.82rem; font-weight:700; color:#115E59;">Energy Captured: ${clean(pct, 2)}% (Pythagoras: $\\|\\hat{f}_N\\|^2 + \\|e\\|^2 = \\|f\\|^2$)</div>
          <div class="m7-energy-bar"><div class="m7-energy-fill" style="width:${pct}%;"></div></div>
          ${taylorErrSq !== null ? `
            <div style="margin-top:6px; font-size:0.82rem; padding:6px 9px; background:#FAF5FF; border:1px solid #E9D5FF; border-radius:7px; color:#581C87;">
              <strong>L² Projection vs. Taylor (deg ${legDeg}):</strong><br/>
              Taylor $L^2$ error $= ${clean(taylorErrSq, 4)}$ vs. Legendre $L^2$ error $= ${clean(errSq, 4)}$
              (${taylorErrSq > errSq + 1e-6 ? `<strong>${clean(taylorErrSq / Math.max(errSq, 1e-9), 1)}× smaller error!</strong>` : 'equal'})
            </div>
          ` : `
            <div style="margin-top:6px; font-size:0.81rem; color:var(--text-muted);">
              Notice: Increasing $N$ <strong>never alters</strong> existing coefficients $c_0, \\dots, c_k$ because $\\{P_k\\}$ is an orthogonal basis!
            </div>
          `}
        `;
        renderMath(statsEl);
      }

      container.querySelectorAll('[data-lf]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-lf]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          activeLf = btn.dataset.lf;
          drawLegLab();
        });
      });
      container.querySelector('#m7Deg').addEventListener('input', e => {
        legDeg = parseInt(e.target.value, 10);
        container.querySelector('#m7DegLbl').textContent = `Degree N = ${legDeg}`;
        drawLegLab();
      });
      container.querySelector('#m7ShowTaylor').addEventListener('change', e => {
        showTaylor = e.target.checked; drawLegLab();
      });
      container.querySelector('#m7ShowBasis').addEventListener('change', e => {
        showBasisOnly = e.target.checked; drawLegLab();
      });

      // =====================================================================
      // LAB 3: FOURIER APPROXIMATION STUDIO IN C[-pi, pi]
      // =====================================================================
      const cFour = container.querySelector('#m7CanvasFour');
      const ctxFour = cFour.getContext('2d');
      const cSpec = container.querySelector('#m7CanvasSpec');
      const ctxSpec = cSpec.getContext('2d');

      const PI = Math.PI;
      const fourSignals = {
        square: {
          label: 'Square Wave',
          fn: t => (t >= 0 ? 1 : -1),
          splits: [0],
          yRange: [-1.5, 1.5],
          note: 'Odd function ($a_k = 0$): only odd sine harmonics $b_k = \\frac{4}{\\pi k}$ survive! Notice the ~9% overshoot near $t=0, \\pm\\pi$ (the <strong>Gibbs phenomenon</strong>).'
        },
        sawtooth: {
          label: 'Sawtooth Wave',
          fn: t => t / PI,
          splits: [],
          yRange: [-1.45, 1.45],
          note: 'All sine harmonics $b_k = \\frac{2(-1)^{k+1}}{\\pi k}$ decay slowly as $O(1/k)$ due to the jump discontinuity at $\\pm\\pi$.'
        },
        triangle: {
          label: 'Triangle Wave',
          fn: t => 1 - (2 * Math.abs(t)) / PI,
          splits: [0],
          yRange: [-1.35, 1.35],
          note: 'Continuous even function ($b_k = 0$): cosine coefficients decay rapidly as $O(1/k^2)$, so just $N=3$ captures $>99.8\\%$ of the energy!'
        },
        rectified: {
          label: 'Half-Wave Rectified Sine',
          fn: t => (t > 0 ? Math.sin(t) : 0),
          splits: [0],
          yRange: [-0.35, 1.35],
          note: 'Classic power-electronics signal: combines a DC offset ($a_0/2 = 1/\\pi$), fundamental sine ($b_1 = 1/2$), and even cosine harmonics!'
        },
        pulse: {
          label: 'Symmetric Pulse',
          fn: t => (Math.abs(t) <= PI / 2 ? 1 : 0),
          splits: [-PI / 2, PI / 2],
          yRange: [-0.35, 1.35],
          note: 'Even function ($b_k = 0$): $a_0/2 = 1/2$ and $a_k = \\frac{2\\sin(k\\pi/2)}{k\\pi}$ alternates signs on odd harmonics.'
        }
      };

      let activeFf = 'square';
      let fourN = 3;
      let animTimer = null;

      function computeFourier(sig, maxN) {
        const a0 = (1 / PI) * integrateSplit(sig.fn, -PI, PI, sig.splits, 140);
        const ak = [a0], bk = [0];
        for (let k = 1; k <= maxN; k++) {
          const a = (1 / PI) * integrateSplit(t => sig.fn(t) * Math.cos(k * t), -PI, PI, sig.splits, 140);
          const b = (1 / PI) * integrateSplit(t => sig.fn(t) * Math.sin(k * t), -PI, PI, sig.splits, 140);
          ak.push(Math.abs(a) < 1e-10 ? 0 : a);
          bk.push(Math.abs(b) < 1e-10 ? 0 : b);
        }
        const normSq = integrateSplit(t => sig.fn(t) * sig.fn(t), -PI, PI, sig.splits, 140);
        return { ak, bk, normSq };
      }

      // Cache up to N=25 per preset so sliding/animating is silky-smooth 60fps
      const fourCache = {};
      function getFourierData(key) {
        if (!fourCache[key]) fourCache[key] = computeFourier(fourSignals[key], 25);
        return fourCache[key];
      }

      function drawFourLab() {
        if (!cFour.isConnected) {
          if (animTimer) { clearInterval(animTimer); animTimer = null; }
          return;
        }
        const sig = fourSignals[activeFf];
        const data = getFourierData(activeFf);

        // --- 1. Main Waveform Canvas ---
        const dpr = window.devicePixelRatio || 1;
        const r = cFour.getBoundingClientRect();
        const w = r.width || 440, h = r.height || 360;
        cFour.width = Math.round(w * dpr);
        cFour.height = Math.round(h * dpr);
        ctxFour.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctxFour.clearRect(0, 0, w, h);

        const padL = 38, padR = 16, padT = 24, padB = 30;
        const plotW = w - padL - padR, plotH = h - padT - padB;
        const [yMin, yMax] = sig.yRange;
        const mapX = t => padL + ((t + PI) / (2 * PI)) * plotW;
        const mapY = y => padT + ((yMax - y) / (yMax - yMin)) * plotH;

        // Grid & axes
        ctxFour.strokeStyle = '#F1F5F9'; ctxFour.lineWidth = 1;
        [-PI, -PI / 2, 0, PI / 2, PI].forEach(t => {
          const x = mapX(t);
          ctxFour.beginPath(); ctxFour.moveTo(x, padT); ctxFour.lineTo(x, h - padB); ctxFour.stroke();
        });
        ctxFour.strokeStyle = '#94A3B8'; ctxFour.lineWidth = 1.3;
        if (yMin <= 0 && yMax >= 0) {
          ctxFour.beginPath(); ctxFour.moveTo(padL, mapY(0)); ctxFour.lineTo(w - padR, mapY(0)); ctxFour.stroke();
        }
        ctxFour.beginPath(); ctxFour.moveTo(mapX(0), padT); ctxFour.lineTo(mapX(0), h - padB); ctxFour.stroke();

        ctxFour.fillStyle = '#64748B'; ctxFour.font = '11px sans-serif';
        [[-PI, '-π'], [-PI / 2, '-π/2'], [0, '0'], [PI / 2, 'π/2'], [PI, 'π']].forEach(([t, lbl]) => {
          ctxFour.fillText(lbl, mapX(t) - 10, h - 10);
        });

        const evalSN = t => {
          let s = data.ak[0] / 2;
          for (let k = 1; k <= fourN; k++) s += data.ak[k] * Math.cos(k * t) + data.bk[k] * Math.sin(k * t);
          return s;
        };

        const samples = 320;
        // Shaded L2 error region
        ctxFour.save();
        ctxFour.fillStyle = 'rgba(217, 119, 6, 0.14)';
        ctxFour.beginPath();
        for (let i = 0; i <= samples; i++) {
          const t = -PI + (2 * PI * i) / samples;
          const x = mapX(t), y = mapY(sig.fn(t));
          if (i === 0) ctxFour.moveTo(x, y); else ctxFour.lineTo(x, y);
        }
        for (let i = samples; i >= 0; i--) {
          const t = -PI + (2 * PI * i) / samples;
          ctxFour.lineTo(mapX(t), mapY(evalSN(t)));
        }
        ctxFour.closePath(); ctxFour.fill();
        ctxFour.restore();

        // Target f(t)
        ctxFour.save();
        ctxFour.strokeStyle = '#1D4ED8'; ctxFour.lineWidth = 2.4; ctxFour.setLineDash([5, 3]);
        ctxFour.beginPath();
        for (let i = 0; i <= samples; i++) {
          const t = -PI + (2 * PI * i) / samples;
          const x = mapX(t), y = mapY(sig.fn(t));
          const jump = (activeFf === 'square' && Math.abs(t) < 0.025) ||
                       (activeFf === 'pulse' && (Math.abs(Math.abs(t) - PI / 2) < 0.025));
          if (i === 0 || jump) ctxFour.moveTo(x, y); else ctxFour.lineTo(x, y);
        }
        ctxFour.stroke();
        ctxFour.restore();

        // Fourier Projection S_N(t)
        ctxFour.save();
        ctxFour.strokeStyle = '#0F766E'; ctxFour.lineWidth = 2.8;
        ctxFour.beginPath();
        for (let i = 0; i <= samples; i++) {
          const t = -PI + (2 * PI * i) / samples;
          const x = mapX(t), y = mapY(evalSN(t));
          if (i === 0) ctxFour.moveTo(x, y); else ctxFour.lineTo(x, y);
        }
        ctxFour.stroke();
        ctxFour.restore();

        ctxFour.fillStyle = '#0F172A'; ctxFour.font = 'bold 11px sans-serif';
        ctxFour.fillText(`Blue dashed: f(t) · Teal: Fourier Projection S_${fourN}(t)`, padL, 15);

        // --- 2. Spectrum Bar Chart Canvas ---
        const rs = cSpec.getBoundingClientRect();
        const sw = rs.width || 400, sh = rs.height || 220;
        cSpec.width = Math.round(sw * dpr);
        cSpec.height = Math.round(sh * dpr);
        ctxSpec.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctxSpec.clearRect(0, 0, sw, sh);

        const maxBars = Math.min(15, Math.max(8, fourN));
        const maxAmp = Math.max(
          0.5,
          Math.abs(data.ak[0] / 2),
          ...data.ak.slice(1, maxBars + 1).map(Math.abs),
          ...data.bk.slice(1, maxBars + 1).map(Math.abs)
        );
        const sPadL = 30, sPadR = 10, sPadT = 16, sPadB = 26;
        const sPlotW = sw - sPadL - sPadR, sPlotH = sh - sPadT - sPadB;
        const slotW = sPlotW / (maxBars + 1);

        ctxSpec.strokeStyle = '#CBD5E1'; ctxSpec.lineWidth = 1.2;
        ctxSpec.beginPath(); ctxSpec.moveTo(sPadL, sh - sPadB); ctxSpec.lineTo(sw - sPadR, sh - sPadB); ctxSpec.stroke();

        for (let k = 0; k <= maxBars; k++) {
          const x0 = sPadL + k * slotW;
          const active = k <= fourN;
          const aVal = k === 0 ? Math.abs(data.ak[0] / 2) : Math.abs(data.ak[k]);
          const bVal = k === 0 ? 0 : Math.abs(data.bk[k]);
          const barW = Math.max(4, slotW * 0.34);

          const hA = (aVal / maxAmp) * sPlotH;
          const hB = (bVal / maxAmp) * sPlotH;
          if (hA > 0.5) {
            ctxSpec.fillStyle = active ? '#0F766E' : '#CBD5E1';
            ctxSpec.fillRect(x0 + 2, sh - sPadB - hA, barW, hA);
          }
          if (hB > 0.5) {
            ctxSpec.fillStyle = active ? '#DC2626' : '#E2E8F0';
            ctxSpec.fillRect(x0 + 2 + barW + 1, sh - sPadB - hB, barW, hB);
          }
          ctxSpec.fillStyle = active ? '#334155' : '#94A3B8';
          ctxSpec.font = '10px sans-serif';
          ctxSpec.fillText(k === 0 ? 'DC' : String(k), x0 + 3, sh - 9);
        }

        // --- 3. Parseval Energy & Formula Readout ---
        // Parseval: ||S_N||^2 = (pi/2)*a0^2 + pi * sum_{k=1}^N (a_k^2 + b_k^2)
        let projEnergy = (PI / 2) * data.ak[0] * data.ak[0];
        for (let k = 1; k <= fourN; k++) projEnergy += PI * (data.ak[k] * data.ak[k] + data.bk[k] * data.bk[k]);
        const errSq = Math.max(0, data.normSq - projEnergy);
        const pct = data.normSq > 1e-12 ? Math.min(100, (projEnergy / data.normSq) * 100) : 100;

        const fStats = container.querySelector('#m7FourStats');
        fStats.innerHTML = `
          <div style="font-size:0.83rem; font-weight:700; color:#115E59;">
            Parseval Energy Captured: ${clean(pct, 2)}% &nbsp;($\\|f - S_${fourN}\\|_{L^2}^2 = ${clean(errSq, 4)}$)
          </div>
          <div class="m7-energy-bar"><div class="m7-energy-fill" style="width:${pct}%;"></div></div>
          <div style="font-size:0.82rem; color:#334155; margin-top:4px;">${sig.note}</div>
        `;
        renderMath(fStats);
      }

      container.querySelectorAll('[data-ff]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-ff]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          activeFf = btn.dataset.ff;
          drawFourLab();
        });
      });
      container.querySelector('#m7Harm').addEventListener('input', e => {
        fourN = parseInt(e.target.value, 10);
        container.querySelector('#m7HarmLbl').textContent = `Harmonics N = ${fourN}`;
        drawFourLab();
      });
      const animBtn = container.querySelector('#m7AnimBtn');
      animBtn.addEventListener('click', () => {
        if (animTimer) {
          clearInterval(animTimer); animTimer = null;
          animBtn.textContent = '▶ Animate N';
        } else {
          fourN = 1;
          container.querySelector('#m7Harm').value = 1;
          container.querySelector('#m7HarmLbl').textContent = 'Harmonics N = 1';
          drawFourLab();
          animBtn.textContent = '⏸ Pause';
          animTimer = setInterval(() => {
            if (!cFour.isConnected) { clearInterval(animTimer); animTimer = null; return; }
            fourN = fourN >= 25 ? 1 : fourN + 1;
            container.querySelector('#m7Harm').value = fourN;
            container.querySelector('#m7HarmLbl').textContent = `Harmonics N = ${fourN}`;
            drawFourLab();
          }, 340);
        }
      });

      // =====================================================================
      // SECTION E: QUICK CONCEPT CHECK
      // =====================================================================
      const explanations = {
        1: 'A nonzero degree-2 polynomial can vanish at two points! Specifically, p(t) = t(t-1) ≠ 0 has p(0) = p(1) = 0, so ⟨p, p⟩ = 0 — violating Axiom 4 (Positivity). For P_n, a point-evaluation inner product requires at least n + 1 distinct sample points!',
        2: 'We need ⟨u, v⟩ = 4(1)(v₁) + 1(2)(v₂) = 4v₁ + 2v₂ = 0, so v₂ = -2v₁. Checking (1, -2)ᵀ gives 4(1) + 2(-2) = 0 ✓! Notice the Euclidean perpendicular (-2, 1)ᵀ fails: 4(-2) + 2(1) = -6 ≠ 0.',
        3: 'Because {P₀, P₁, P₂, P₃} is an orthogonal basis, each coordinate c_k = ⟨f, P_k⟩ / ⟨P_k, P_k⟩ depends ONLY on f and P_k — completely independent of the other basis polynomials! That is why orthogonal polynomials and Fourier series are so prized in numerical computing.',
        4: 'Set g(t) = 1. Then ⟨f, 1⟩ = ∫₀¹ f(t) dt, ‖1‖² = ∫₀¹ 1² dt = 1, and ‖f‖² = ∫₀¹ f(t)² dt. Squaring Cauchy–Schwarz |⟨f, 1⟩| ≤ ‖f‖‖1‖ immediately yields (∫₀¹ f(t) dt)² ≤ ∫₀¹ f(t)² dt!'
      };

      container.querySelectorAll('.m7-quiz-item').forEach(item => {
        const qNum = item.dataset.q;
        const ans = item.dataset.ans;
        const fb = container.querySelector(`#m7Fb${qNum}`);
        item.querySelectorAll('.m7-opt-btn').forEach(btn => {
          btn.addEventListener('click', () => {
            item.querySelectorAll('.m7-opt-btn').forEach(b => b.classList.remove('correct', 'wrong'));
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
        drawWLab(); drawLegLab(); drawFourLab();
      });

      // Initial renders
      syncWSliders();
      drawWLab();
      drawLegLab();
      drawFourLab();
    }
  });
})();
