/**
 * AMA1751 Linear Algebra — Topic 6: Inner Product
 * Module 10: Capstone Synthesis — Orthogonal Matrices (Isometries),
 *            Continuous Inner Product Spaces L^2[-1,1] (Legendre Projection),
 *            and Topic 6 Comprehensive Timed Exam Simulator
 */
(function () {
  function clean(x, d = 3) {
    if (!Number.isFinite(x)) return 'NaN';
    if (Math.abs(x) < 1e-10) return '0';
    return String(Number(x.toFixed(d)));
  }

  function dot2(a, b) { return a[0] * b[0] + a[1] * b[1]; }
  function norm2(a) { return Math.hypot(a[0], a[1]); }
  function matVec2(M, v) {
    return [M[0][0] * v[0] + M[0][1] * v[1], M[1][0] * v[0] + M[1][1] * v[1]];
  }
  function angleDeg2(a, b) {
    const na = norm2(a), nb = norm2(b);
    if (na < 1e-9 || nb < 1e-9) return 0;
    const c = Math.max(-1, Math.min(1, dot2(a, b) / (na * nb)));
    return (Math.acos(c) * 180) / Math.PI;
  }

  // =========================================================================
  // LAB 2 MATH: LEGENDRE POLYNOMIALS & COMPOSITE SIMPSON ON [-1, 1]
  // =========================================================================
  function legendreP(n, x) {
    if (n === 0) return 1;
    if (n === 1) return x;
    let p0 = 1, p1 = x;
    for (let k = 1; k < n; k++) {
      const pk = ((2 * k + 1) * x * p1 - k * p0) / (k + 1);
      p0 = p1;
      p1 = pk;
    }
    return p1;
  }

  // Composite Simpson on [-1, 1] with M = 400 (even, so x = 0 is a node!)
  function innerL2(f, g) {
    const M = 400;
    const h = 2 / M;
    let sum = f(-1) * g(-1) + f(1) * g(1);
    for (let i = 1; i < M; i++) {
      const x = -1 + i * h;
      sum += (i % 2 === 0 ? 2 : 4) * f(x) * g(x);
    }
    return (h / 3) * sum;
  }

  const L2_FUNCS = {
    runge: {
      name: 'Runge Function',
      tex: 'f(x) = \\dfrac{1}{1 + 4x^2}',
      fn: x => 1 / (1 + 4 * x * x),
      hasTaylor: true,
      taylor: (n, x) => {
        let s = 0;
        for (let m = 0; 2 * m <= n; m++) s += Math.pow(-4, m) * Math.pow(x, 2 * m);
        return s;
      },
      note: 'Its Taylor series diverges for $|x| \\ge 0.5$! Meanwhile, the $L^2$-orthogonal projection onto $\\mathcal{P}_n$ converges cleanly across the entire interval $[-1,1]$.'
    },
    abs: {
      name: 'Absolute Value',
      tex: 'f(x) = |x|',
      fn: x => Math.abs(x),
      hasTaylor: false,
      taylor: null,
      note: '$|x|$ is not differentiable at $x=0$, so no Maclaurin polynomial exists—yet the inner product $\\langle f, P_k \\rangle$ still delivers the unique best $L^2$ polynomial approximation!'
    },
    sinpi: {
      name: 'Sine Wave',
      tex: 'f(x) = \\sin(\\pi x)',
      fn: x => Math.sin(Math.PI * x),
      hasTaylor: true,
      taylor: (n, x) => {
        let s = 0, fact = 1;
        for (let k = 1; k <= n; k++) {
          fact *= k;
          if (k % 2 === 1) {
            const m = (k - 1) / 2;
            s += ((m % 2 === 0 ? 1 : -1) * Math.pow(Math.PI * x, k)) / fact;
          }
        }
        return s;
      },
      note: 'Taylor matches $\\sin(\\pi x)$ tightly near $x=0$ but peels off dramatically near $x = \\pm 1$; the Legendre projection spreads error optimally over all of $[-1,1]$.'
    },
    exp: {
      name: 'Exponential',
      tex: 'f(x) = e^x',
      fn: x => Math.exp(x),
      hasTaylor: true,
      taylor: (n, x) => {
        let s = 1, term = 1;
        for (let k = 1; k <= n; k++) { term *= x / k; s += term; }
        return s;
      },
      note: 'Compare degree $n=1$ or $n=2$: Taylor uses the tangent line $1+x$, whereas the $L^2$ projection shifts up slightly to minimize global squared error $\\int_{-1}^1 (e^x - p(x))^2 dx$.'
    },
    step: {
      name: 'Sign Step',
      tex: 'f(x) = \\operatorname{sgn}(x)',
      fn: x => (x < 0 ? -1 : x > 0 ? 1 : 0),
      hasTaylor: false,
      taylor: null,
      note: 'Even for a discontinuous step function, orthogonal projection in $L^2[-1,1]$ exists uniquely (notice the odd-only Legendre coefficients and Gibbs overshoot at $x=0$!).'
    }
  };

  // =========================================================================
  // SECTION C: TOPIC 6 COMPREHENSIVE EXAM QUESTION BANK (10 VERIFIED Qs)
  // =========================================================================
  const EXAM_BANK = [
    {
      id: 1,
      pillar: 'Inner Products & Angles',
      q: 'Let $\\mathbf{u} = (1, -2, 2)^T$ and $\\mathbf{v} = (2, 1, -2)^T$ in $\\mathbb{R}^3$. What is $\\cos\\theta$ between $\\mathbf{u}$ and $\\mathbf{v}$?',
      opts: ['$-\\dfrac{4}{9}$', '$\\dfrac{4}{9}$', '$-\\dfrac{2}{3}$', '$0$'],
      ans: 0,
      exp: '$\\mathbf{u}\\cdot\\mathbf{v} = (1)(2)+(-2)(1)+(2)(-2) = -4$, $\\|\\mathbf{u}\\| = \\sqrt{1+4+4}=3$, and $\\|\\mathbf{v}\\| = \\sqrt{4+1+4}=3$. Hence $\\cos\\theta = \\dfrac{-4}{3\\cdot 3} = -\\dfrac{4}{9}$.'
    },
    {
      id: 2,
      pillar: 'Projections & Orthonormal Bases',
      q: 'Decompose $\\mathbf{y} = (3, 4)^T$ as $\\mathbf{y} = \\hat{\\mathbf{y}} + \\mathbf{z}$ where $\\hat{\\mathbf{y}}$ is along $\\mathbf{u} = (2, 1)^T$ and $\\mathbf{z} \\perp \\mathbf{u}$. What is $\\mathbf{z}$?',
      opts: ['$(4, 2)^T$', '$(-1, 2)^T$', '$(1, -2)^T$', '$(-2, 4)^T$'],
      ans: 1,
      exp: '$\\hat{\\mathbf{y}} = \\dfrac{\\mathbf{y}\\cdot\\mathbf{u}}{\\mathbf{u}\\cdot\\mathbf{u}}\\mathbf{u} = \\dfrac{10}{5}(2,1)^T = (4,2)^T$. Thus $\\mathbf{z} = \\mathbf{y} - \\hat{\\mathbf{y}} = (3-4, 4-2)^T = (-1, 2)^T$ (check: $\\mathbf{z}\\cdot\\mathbf{u} = -2+2 = 0$).'
    },
    {
      id: 3,
      pillar: 'Projections & Orthonormal Bases',
      q: 'Let $W = \\operatorname{span}\\{\\mathbf{q}_1, \\mathbf{q}_2\\} \\subset \\mathbb{R}^3$ where $\\mathbf{q}_1 = \\frac{1}{3}(2, 2, -1)^T$ and $\\mathbf{q}_2 = \\frac{1}{3}(2, -1, 2)^T$ are orthonormal. What is the shortest distance from $\\mathbf{y} = (3, 0, 3)^T$ to $W$?',
      opts: ['$0$', '$1$', '$\\sqrt{2}$', '$3$'],
      ans: 1,
      exp: 'Compute orthonormal coordinates: $c_1 = \\mathbf{y}\\cdot\\mathbf{q}_1 = \\frac{6-3}{3} = 1$, $c_2 = \\mathbf{y}\\cdot\\mathbf{q}_2 = \\frac{6+6}{3} = 4$. By Pythagoras, $\\|\\mathbf{y} - \\hat{\\mathbf{y}}\\|^2 = \\|\\mathbf{y}\\|^2 - (c_1^2 + c_2^2) = 18 - (1 + 16) = 1$, so the distance is $\\sqrt{1} = 1$.'
    },
    {
      id: 4,
      pillar: 'Gram–Schmidt, QR & Stability',
      q: 'Compute the upper-triangular factor $R$ (with positive diagonal entries) in the $QR$ factorization of $A = \\begin{bmatrix} 3 & 5 \\\\ 4 & 0 \\end{bmatrix}$.',
      opts: [
        '$\\begin{bmatrix} 5 & 3 \\\\ 0 & 4 \\end{bmatrix}$',
        '$\\begin{bmatrix} 5 & 4 \\\\ 0 & 3 \\end{bmatrix}$',
        '$\\begin{bmatrix} 3 & 5 \\\\ 0 & 4 \\end{bmatrix}$',
        '$\\begin{bmatrix} 5 & 3 \\\\ 0 & -4 \\end{bmatrix}$'
      ],
      ans: 0,
      exp: '$r_{11} = \\|(3,4)^T\\| = 5$, $\\mathbf{q}_1 = \\frac{1}{5}(3,4)^T$. Then $r_{12} = \\mathbf{q}_1\\cdot(5,0)^T = 3$, $\\mathbf{v}_2 = (5,0)^T - 3\\mathbf{q}_1 = \\frac{4}{5}(4,-3)^T$, so $r_{22} = \\|\\mathbf{v}_2\\| = 4$.'
    },
    {
      id: 5,
      pillar: 'Subspaces & Least Squares',
      q: 'Suppose $A$ is a $5 \\times 7$ matrix with $\\operatorname{rank}(A) = 3$. What are $\\dim(\\operatorname{Col}(A)^\\perp)$ and $\\dim(\\operatorname{Row}(A)^\\perp)$?',
      opts: [
        '$2$ and $4$',
        '$4$ and $2$',
        '$3$ and $4$',
        '$2$ and $3$'
      ],
      ans: 0,
      exp: '$\\operatorname{Col}(A) \\subseteq \\mathbb{R}^5$ has dimension $3$, so $\\dim(\\operatorname{Col}(A)^\\perp) = \\dim(\\operatorname{Null}(A^T)) = 5 - 3 = 2$. Similarly, $\\operatorname{Row}(A) \\subseteq \\mathbb{R}^7$ has dimension $3$, so $\\dim(\\operatorname{Row}(A)^\\perp) = \\dim(\\operatorname{Null}(A)) = 7 - 3 = 4$.'
    },
    {
      id: 6,
      pillar: 'Subspaces & Least Squares',
      q: 'Find the least-squares solution $\\hat{\\mathbf{x}}$ to $A\\mathbf{x} = \\mathbf{b}$ where $A = \\begin{bmatrix} 1 & 0 \\\\ 1 & 1 \\\\ 1 & 2 \\end{bmatrix}$ and $\\mathbf{b} = \\begin{bmatrix} 2 \\\\ 1 \\\\ 6 \\end{bmatrix}$.',
      opts: [
        '$\\hat{\\mathbf{x}} = (1, 2)^T$',
        '$\\hat{\\mathbf{x}} = (2, 1)^T$',
        '$\\hat{\\mathbf{x}} = (0, 3)^T$',
        '$\\hat{\\mathbf{x}} = (1, 3)^T$'
      ],
      ans: 0,
      exp: 'Normal equations: $A^T A = \\begin{bmatrix} 3 & 3 \\\\ 3 & 5 \\end{bmatrix}$ and $A^T\\mathbf{b} = \\begin{bmatrix} 9 \\\\ 13 \\end{bmatrix}$. Subtracting row 1 from row 2 gives $2\\hat{x}_2 = 4 \\Rightarrow \\hat{x}_2 = 2$, and $3\\hat{x}_1 + 6 = 9 \\Rightarrow \\hat{x}_1 = 1$.'
    },
    {
      id: 7,
      pillar: 'Orthogonal Matrices & Function Spaces',
      q: 'Which statement is <strong>always true</strong> for every $n \\times n$ orthogonal matrix $Q$ (satisfying $Q^T Q = I_n$)?',
      opts: [
        'Every eigenvalue of $Q$ is $+1$ or $-1$.',
        '$Q$ preserves lengths and inner products, $Q^{-1} = Q^T$, and $\\det(Q) = \\pm 1$.',
        'Every square matrix with $\\det(A) = \\pm 1$ is orthogonal.',
        '$Q$ must be symmetric ($Q^T = Q$).'
      ],
      ans: 1,
      exp: 'An orthogonal matrix satisfies $(Q\\mathbf{u})\\cdot(Q\\mathbf{v}) = \\mathbf{u}^T Q^T Q \\mathbf{v} = \\mathbf{u}\\cdot\\mathbf{v}$, $Q^{-1}=Q^T$, and $(\\det Q)^2 = \\det(Q^T Q) = 1$. Note that shears have $\\det = 1$ without being orthogonal, and rotations are generally not symmetric and can have complex eigenvalues $e^{\\pm i\\theta}$.'
    },
    {
      id: 8,
      pillar: 'Orthogonal Matrices & Function Spaces',
      q: 'In the inner product space $C[-1, 1]$ with $\\langle f, g \\rangle = \\int_{-1}^1 f(x)g(x)\\,dx$, for which value of $c$ is $p(x) = x^2 - c$ orthogonal to $q(x) = 1$?',
      opts: [
        '$c = 0$',
        '$c = \\dfrac{1}{2}$',
        '$c = \\dfrac{1}{3}$',
        '$c = \\dfrac{2}{3}$'
      ],
      ans: 2,
      exp: '$\\langle x^2 - c, 1 \\rangle = \\int_{-1}^1 (x^2 - c)\\,dx = \\left[\\frac{x^3}{3} - cx\\right]_{-1}^1 = \\frac{2}{3} - 2c = 0 \\implies c = \\frac{1}{3}$. (Notice $x^2 - \\frac{1}{3} = \\frac{2}{3} P_2(x)$, the degree-2 Legendre polynomial!)'
    },
    {
      id: 9,
      pillar: 'Inner Products & Angles',
      q: 'Let $P = Q Q^T$, where $Q$ is an $m \\times k$ matrix with orthonormal columns and $1 \\le k < m$. Which statement is <strong>FALSE</strong>?',
      opts: [
        '$P^2 = P$ (idempotent)',
        '$P^T = P$ (symmetric)',
        '$P$ is invertible',
        '$\\|P\\mathbf{y}\\| \\le \\|\\mathbf{y}\\|$ for all $\\mathbf{y} \\in \\mathbb{R}^m$'
      ],
      ans: 2,
      exp: 'Since $\\operatorname{rank}(P) = k < m$, $P$ is <strong>singular</strong> (not invertible): every nonzero vector $\\mathbf{z} \\in \\operatorname{Col}(Q)^\\perp$ satisfies $P\\mathbf{z} = \\mathbf{0}$. The other three properties hold for every orthogonal projection matrix.'
    },
    {
      id: 10,
      pillar: 'Gram–Schmidt, QR & Stability',
      q: 'Why do production numerical libraries (like NumPy and MATLAB) prefer <strong>Householder QR</strong> ($R\\hat{\\mathbf{x}} = Q^T\\mathbf{b}$) over the normal equations ($A^T A\\hat{\\mathbf{x}} = A^T\\mathbf{b}$) for ill-conditioned least-squares problems?',
      opts: [
        'Because normal equations only work when $A$ is square.',
        'Because forming $A^T A$ squares the condition number ($\\kappa(A^T A) = \\kappa(A)^2$), risking catastrophic floating-point digit loss, whereas Householder QR works directly on $A$.',
        'Because Householder reflections always make $R$ diagonal.',
        'Because $A^T A$ is never symmetric.'
      ],
      ans: 1,
      exp: 'Squaring the condition number can wipe out all 16 digits of IEEE-754 double precision as soon as $\\kappa(A) \\approx 10^8$. Householder reflections are orthogonal transformations that preserve norms and avoid forming $A^T A$.'
    }
  ];

  // =========================================================================
  // MODULE 10 REGISTRATION
  // =========================================================================
  AMA1751.registerModule(10, {
    render(meta) {
      const slideLabel = (meta && meta.slides) ? meta.slides : 'Capstone Synthesis · Topic 6 Mastery';
      return `
        <style>
          .m10-controls {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            margin-bottom: 12px;
            align-items: center;
          }
          .m10-chip {
            padding: 6px 12px;
            border-radius: 999px;
            border: 1px solid var(--border);
            background: #F8FAFC;
            font-size: 0.82rem;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.15s;
          }
          .m10-chip:hover, .m10-chip.active {
            background: var(--polyu-red-soft);
            border-color: #E8B4C4;
            color: var(--polyu-red-dark);
          }
          .m10-btn {
            padding: 7px 13px;
            border-radius: 8px;
            border: 1px solid #CBD5E1;
            background: #FFFFFF;
            font-weight: 700;
            font-size: 0.83rem;
            cursor: pointer;
          }
          .m10-btn:hover { background: #F1F5F9; }
          .m10-btn.primary {
            background: var(--polyu-red-dark);
            border-color: var(--polyu-red-dark);
            color: #FFFFFF;
          }
          .m10-btn.primary:hover { opacity: 0.92; }
          .m10-panel {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 14px 16px;
            font-size: 0.89rem;
            line-height: 1.55;
          }
          .m10-canvas-box {
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 8px;
          }
          .m10-canvas {
            display: block;
            width: 100%;
          }
          .m10-slider-row {
            display: flex;
            align-items: center;
            gap: 10px;
            margin: 6px 0;
            font-size: 0.85rem;
          }
          .m10-slider-row label {
            min-width: 140px;
            font-weight: 700;
            color: #334155;
          }
          .m10-slider-row input[type="range"] {
            flex: 1;
            accent-color: #0F766E;
          }
          .m10-slider-val {
            min-width: 72px;
            text-align: right;
            font-weight: 800;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
            color: #0F172A;
          }
          .m10-stat-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(128px, 1fr));
            gap: 8px;
            margin: 10px 0;
          }
          .m10-stat-card {
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 8px 10px;
            text-align: center;
          }
          .m10-stat-lbl {
            font-size: 0.69rem;
            font-weight: 700;
            text-transform: uppercase;
            color: var(--text-muted);
          }
          .m10-stat-val {
            font-size: 0.94rem;
            font-weight: 800;
            margin-top: 2px;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
            color: #0F172A;
          }
          .m10-pill {
            display: inline-block;
            padding: 2px 8px;
            border-radius: 999px;
            font-size: 0.74rem;
            font-weight: 800;
          }
          .m10-pill.ok { background: #DCFCE7; color: #166534; }
          .m10-pill.bad { background: #FEE2E2; color: #991B1B; }
          .m10-pill.info { background: #DBEAFE; color: #1E40AF; }
          /* Exam Simulator Styles */
          .m10-exam-bar {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            align-items: center;
            justify-content: space-between;
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 10px 14px;
            margin-bottom: 12px;
          }
          .m10-timer {
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
            font-weight: 800;
            font-size: 1rem;
            color: #0F172A;
            background: #FFFFFF;
            border: 1px solid #CBD5E1;
            padding: 4px 10px;
            border-radius: 8px;
          }
          .m10-qnav {
            display: flex;
            flex-wrap: wrap;
            gap: 6px;
            margin-bottom: 12px;
          }
          .m10-qdot {
            width: 34px;
            height: 30px;
            border-radius: 7px;
            border: 1px solid #CBD5E1;
            background: #FFFFFF;
            font-weight: 700;
            font-size: 0.8rem;
            cursor: pointer;
          }
          .m10-qdot.active { border-color: #0F172A; box-shadow: 0 0 0 2px rgba(15,23,42,0.18); }
          .m10-qdot.answered { background: #E0F2FE; border-color: #0284C7; color: #0369A1; }
          .m10-qdot.correct { background: #DCFCE7; border-color: #16A34A; color: #166534; }
          .m10-qdot.wrong { background: #FEE2E2; border-color: #DC2626; color: #991B1B; }
          .m10-opt-btn {
            display: block;
            width: 100%;
            text-align: left;
            padding: 10px 12px;
            margin-top: 8px;
            border: 1px solid #CBD5E1;
            border-radius: 8px;
            background: #FFFFFF;
            cursor: pointer;
            font-size: 0.89rem;
            transition: all 0.12s;
          }
          .m10-opt-btn:hover { background: #F1F5F9; }
          .m10-opt-btn.selected { border-color: #0284C7; background: #E0F2FE; font-weight: 700; }
          .m10-opt-btn.correct { border-color: #16A34A; background: #DCFCE7; color: #166534; font-weight: 700; }
          .m10-opt-btn.wrong { border-color: #DC2626; background: #FEE2E2; color: #991B1B; }
          .m10-bar-track {
            width: 100%;
            height: 9px;
            background: #E2E8F0;
            border-radius: 999px;
            overflow: hidden;
            margin-top: 4px;
          }
          .m10-bar-fill {
            height: 100%;
            border-radius: 999px;
            transition: width 0.25s;
          }
        </style>

        <!-- HEADER -->
        <div class="module-header">
          <span class="slide-ref">Part 10 · ${slideLabel}</span>
          <h2>10. Capstone: Orthogonal Isometries, Function Spaces &amp; Topic 6 Exam Simulator</h2>
          <p>Bring all of <strong>Topic 6 (Inner Product)</strong> together: explore orthogonal transformations as rigid isometries, lift orthogonal projections into continuous function spaces $L^2[-1,1]$, and test your exam readiness across all 5 core pillars.</p>
        </div>

        <!-- =============================================================== -->
        <!-- LAB 1: ORTHOGONAL MATRICES & GEOMETRIC ISOMETRIES O(2)          -->
        <!-- =============================================================== -->
        <div class="card">
          <h3>Interactive Lab 1: Orthogonal Matrices ($Q^T Q = I$) vs. Area-Preserving Shear</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            A square matrix $Q$ is <strong>orthogonal</strong> iff $Q^T Q = I$, which means $Q$ preserves every inner product, length, and angle: $(Q\mathbf{u})\cdot(Q\mathbf{v}) = \mathbf{u}\cdot\mathbf{v}$.
            <strong>Drag $\mathbf{u}$ (teal) or $\mathbf{v}$ (amber)</strong> directly on the canvas, and compare a Rotation, a Reflection, and a Non-Orthogonal Shear!
          </p>

          <div class="m10-controls">
            <strong style="font-size:0.82rem;">Transformation $A$:</strong>
            <button class="m10-chip active" data-iso="rot">Rotation $R(\\theta)$ ($\\det=+1$)</button>
            <button class="m10-chip" data-iso="ref">Reflection $H(\\phi)$ ($\\det=-1$)</button>
            <button class="m10-chip" data-iso="shear">Non-Orthogonal Shear $S(k)$ ($\\det=+1$ Trap!)</button>
          </div>

          <div class="m10-panel" style="margin-bottom:12px;">
            <div class="m10-slider-row">
              <label for="m10IsoSlider" id="m10IsoLabel">Rotation Angle $\\theta$:</label>
              <input type="range" id="m10IsoSlider" min="-180" max="180" step="1" value="50" />
              <span class="m10-slider-val" id="m10IsoVal">50°</span>
            </div>
          </div>

          <div class="grid-2" style="align-items:start;">
            <div class="m10-canvas-box">
              <canvas id="m10IsoCanvas" class="m10-canvas" style="height:300px; cursor:grab;"></canvas>
              <div style="font-size:0.77rem; color:var(--text-muted); text-align:center; margin-top:4px;">
                Drag tips of <strong>u</strong> or <strong>v</strong> · Dashed circle = unit circle · Solid purple curve = image of unit circle under $A$
              </div>
            </div>
            <div class="m10-panel">
              <div id="m10IsoMatrix" style="overflow-x:auto; margin-bottom:8px;"></div>
              <div class="m10-stat-grid" id="m10IsoStats"></div>
              <div id="m10IsoVerdict" style="font-size:0.84rem; margin-top:8px;"></div>
            </div>
          </div>
        </div>

        <!-- =============================================================== -->
        <!-- LAB 2: CONTINUOUS INNER PRODUCT SPACE L^2[-1,1] & LEGENDRE      -->
        <!-- =============================================================== -->
        <div class="card">
          <h3>Interactive Lab 2: Function Inner Product Space $L^2[-1,1]$ &amp; Legendre Orthogonal Projection</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            Inner products extend beyond $\\mathbb{R}^n$! On $C[-1,1]$ with $\\langle f, g \\rangle = \\int_{-1}^1 f(x)g(x)\\,dx$, applying Gram–Schmidt to $\\{1, x, x^2, \\dots\\}$ produces the <strong>Legendre orthogonal polynomials</strong> $\\{P_k(x)\\}$ with $\\|P_k\\|_{L^2}^2 = \\frac{2}{2k+1}$.
            The best polynomial approximation of degree $\\le n$ is simply the orthogonal projection $\\hat{f}_n(x) = \\sum_{k=0}^n \\frac{\\langle f, P_k \\rangle}{\\|P_k\\|^2} P_k(x)$!
          </p>

          <div class="m10-controls">
            <strong style="font-size:0.82rem;">Target $f(x)$:</strong>
            <button class="m10-chip active" data-l2="runge">Runge $\\frac{1}{1+4x^2}$</button>
            <button class="m10-chip" data-l2="abs">$|x|$ (No Taylor!)</button>
            <button class="m10-chip" data-l2="sinpi">$\\sin(\\pi x)$</button>
            <button class="m10-chip" data-l2="exp">$e^x$</button>
            <button class="m10-chip" data-l2="step">Step $\\operatorname{sgn}(x)$</button>
            <label style="margin-left:auto; font-size:0.82rem; font-weight:700; display:flex; align-items:center; gap:5px; cursor:pointer;">
              <input type="checkbox" id="m10ShowTaylor" checked /> Overlay Taylor $T_n(x)$
            </label>
          </div>

          <div class="m10-panel" style="margin-bottom:12px;">
            <div class="m10-slider-row">
              <label for="m10DegSlider">Subspace Degree $\\mathcal{P}_n$:</label>
              <input type="range" id="m10DegSlider" min="0" max="7" step="1" value="4" />
              <span class="m10-slider-val" id="m10DegVal">n = 4</span>
            </div>
          </div>

          <div class="m10-canvas-box" style="margin-bottom:10px;">
            <canvas id="m10L2Canvas" class="m10-canvas" style="height:280px;"></canvas>
          </div>

          <div class="m10-stat-grid" id="m10L2Stats"></div>
          <div class="m10-panel" id="m10L2Readout" style="font-size:0.85rem; margin-bottom:10px;"></div>

          <details class="code-details">
            <summary>Python (NumPy / SciPy) Companion Code — Orthogonal Matrices &amp; Legendre $L^2$ Projection</summary>
            <pre><code>import numpy as np
from scipy.integrate import quad
from numpy.polynomial.legendre import Legendre

# 1) Verify Rotation & Reflection preserve inner products, whereas Shear does not
theta = np.deg2rad(50)
R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
S = np.array([[1.0, 0.8], [0.0, 1.0]])  # det(S) = 1, yet NOT orthogonal!
u, v = np.array([2.0, 0.5]), np.array([0.5, 1.8])
print("u.v =", np.dot(u, v), "| (Ru).(Rv) =", np.dot(R @ u, R @ v), "| (Su).(Sv) =", np.dot(S @ u, S @ v))

# 2) L^2[-1, 1] Orthogonal Projection of f(x) = 1 / (1 + 4x^2) onto P_4 via Legendre Basis
f = lambda x: 1.0 / (1.0 + 4.0 * x**2)
n = 4
coeffs = []
for k in range(n + 1):
    Pk = Legendre.basis(k)
    inner_f_Pk, _ = quad(lambda x: f(x) * Pk(x), -1, 1)
    norm_sq = 2.0 / (2 * k + 1)
    coeffs.append(inner_f_Pk / norm_sq)

f_hat = Legendre(coeffs)  # Best-fit polynomial in L^2[-1, 1]
err_sq, _ = quad(lambda x: (f(x) - f_hat(x))**2, -1, 1)
print("Legendre coefficients a_0..a_4:", np.round(coeffs, 4))
print("L^2 projection error ||f - f_hat||_L2 =", np.sqrt(err_sq))</code></pre>
          </details>
        </div>

        <!-- =============================================================== -->
        <!-- SECTION C: TOPIC 6 COMPREHENSIVE EXAM SIMULATOR (10 QUESTIONS)  -->
        <!-- =============================================================== -->
        <div class="card">
          <h3>Topic 6 Comprehensive Exam Simulator (10 Questions · Mastery Diagnostic)</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            Test yourself across all 5 core Topic 6 pillars. Use <strong>Practice Mode</strong> for instant step-by-step solutions after every click, or switch to <strong>Timed Exam Mode</strong> to simulate a real quiz and grade all 10 questions at once!
          </p>

          <div class="m10-exam-bar">
            <div style="display:flex; flex-wrap:wrap; gap:6px; align-items:center;">
              <button class="m10-chip active" data-emode="practice">Practice Mode</button>
              <button class="m10-chip" data-emode="exam">Timed Exam Mode</button>
              <button class="m10-btn" id="m10ShuffleBtn">⇄ Shuffle</button>
              <button class="m10-btn" id="m10ResetExamBtn">↺ Reset</button>
            </div>
            <div style="display:flex; gap:8px; align-items:center;">
              <span class="m10-timer" id="m10TimerDisplay">00:00</span>
              <button class="m10-btn" id="m10TimerToggleBtn">▶ Start Timer</button>
            </div>
          </div>

          <div class="m10-qnav" id="m10QNav"></div>

          <div class="m10-panel" id="m10QuestionBox" style="background:#FFFFFF;"></div>

          <div style="display:flex; justify-content:space-between; align-items:center; margin-top:10px; gap:8px; flex-wrap:wrap;">
            <div style="display:flex; gap:8px;">
              <button class="m10-btn" id="m10PrevQ">← Previous</button>
              <button class="m10-btn" id="m10NextQ">Next →</button>
            </div>
            <button class="m10-btn primary" id="m10SubmitExamBtn">Grade Full Exam &amp; Show Mastery Report</button>
          </div>

          <div id="m10ReportBox" style="display:none; margin-top:14px;"></div>
        </div>
      `;
    },

    mount(container, { renderMath }) {
      // =====================================================================
      // LAB 1: ORTHOGONAL MATRICES & GEOMETRIC ISOMETRIES O(2)
      // =====================================================================
      let isoMode = 'rot'; // 'rot' | 'ref' | 'shear'
      let uVec = [2.0, 0.5];
      let vVec = [0.6, 1.8];
      let dragging = null;

      const isoSlider = container.querySelector('#m10IsoSlider');
      const isoLabel = container.querySelector('#m10IsoLabel');
      const isoVal = container.querySelector('#m10IsoVal');
      const isoCanvas = container.querySelector('#m10IsoCanvas');
      const isoCtx = isoCanvas.getContext('2d');

      function getIsoMatrix() {
        const val = parseFloat(isoSlider.value);
        if (isoMode === 'rot') {
          const rad = (val * Math.PI) / 180;
          const c = Math.cos(rad), s = Math.sin(rad);
          return { M: [[c, -s], [s, c]], param: val, rad };
        }
        if (isoMode === 'ref') {
          const rad = (val * Math.PI) / 180;
          const c2 = Math.cos(2 * rad), s2 = Math.sin(2 * rad);
          return { M: [[c2, s2], [s2, -c2]], param: val, rad };
        }
        // Shear S(k)
        return { M: [[1, val], [0, 1]], param: val, rad: 0 };
      }

      function getIsoTransformCoords(w, h) {
        const cx = w / 2, cy = h / 2;
        const scale = Math.min(w, h) / 7.4; // fits [-3.7, 3.7]
        return {
          toPx: p => [cx + p[0] * scale, cy - p[1] * scale],
          toMath: (px, py) => [(px - cx) / scale, (cy - py) / scale],
          cx, cy, scale
        };
      }

      function updateIsoLab() {
        const { M, param, rad } = getIsoMatrix();
        isoVal.textContent = isoMode === 'shear' ? `k = ${param.toFixed(2)}` : `${param}°`;

        const Au = matVec2(M, uVec);
        const Av = matVec2(M, vVec);

        const nu = norm2(uVec), nAu = norm2(Au);
        const nv = norm2(vVec), nAv = norm2(Av);
        const dBefore = dot2(uVec, vVec), dAfter = dot2(Au, Av);
        const angBefore = angleDeg2(uVec, vVec), angAfter = angleDeg2(Au, Av);

        // Compute ||M^T M - I||_F and det(M)
        const mtm00 = M[0][0] * M[0][0] + M[1][0] * M[1][0] - 1;
        const mtm01 = M[0][0] * M[0][1] + M[1][0] * M[1][1];
        const mtm11 = M[0][1] * M[0][1] + M[1][1] * M[1][1] - 1;
        const orthErr = Math.sqrt(mtm00 * mtm00 + 2 * mtm01 * mtm01 + mtm11 * mtm11);
        const detM = M[0][0] * M[1][1] - M[0][1] * M[1][0];
        const isOrth = orthErr < 1e-6;

        const matEl = container.querySelector('#m10IsoMatrix');
        matEl.innerHTML = `$$A = \\begin{bmatrix} ${clean(M[0][0], 3)} & ${clean(M[0][1], 3)} \\\\ ${clean(M[1][0], 3)} & ${clean(M[1][1], 3)} \\end{bmatrix}, \\quad \\det(A) = ${clean(detM, 2)}$$`;
        renderMath(matEl);

        const statsEl = container.querySelector('#m10IsoStats');
        statsEl.innerHTML = `
          <div class="m10-stat-card">
            <div class="m10-stat-lbl">$\\|\\mathbf{u}\\| \\to \\|A\\mathbf{u}\\|$</div>
            <div class="m10-stat-val" style="color:${Math.abs(nu - nAu) < 1e-4 ? '#0F766E' : '#DC2626'};">${clean(nu, 2)} → ${clean(nAu, 2)}</div>
          </div>
          <div class="m10-stat-card">
            <div class="m10-stat-lbl">$\\mathbf{u}\\cdot\\mathbf{v} \\to (A\\mathbf{u})\\cdot(A\\mathbf{v})$</div>
            <div class="m10-stat-val" style="color:${Math.abs(dBefore - dAfter) < 1e-4 ? '#0F766E' : '#DC2626'};">${clean(dBefore, 2)} → ${clean(dAfter, 2)}</div>
          </div>
          <div class="m10-stat-card">
            <div class="m10-stat-lbl">$\\theta \\to \\theta_{\\text{after}}$</div>
            <div class="m10-stat-val" style="color:${Math.abs(angBefore - angAfter) < 0.05 ? '#0F766E' : '#DC2626'};">${angBefore.toFixed(1)}° → ${angAfter.toFixed(1)}°</div>
          </div>
          <div class="m10-stat-card">
            <div class="m10-stat-lbl">$\\|A^T A - I\\|_F$</div>
            <div class="m10-stat-val" style="color:${isOrth ? '#0F766E' : '#DC2626'};">${clean(orthErr, 3)}</div>
          </div>
        `;
        renderMath(statsEl);

        const verdictEl = container.querySelector('#m10IsoVerdict');
        if (isoMode === 'rot') {
          verdictEl.innerHTML = `<span class="m10-pill ok">Orthogonal Isometry (Proper Rotation)</span> $Q^T Q = I$ and $\\det(Q) = +1$: preserves lengths, angles, inner products, and orientation! Notice the unit circle maps onto itself.`;
        } else if (isoMode === 'ref') {
          verdictEl.innerHTML = `<span class="m10-pill ok">Orthogonal Isometry (Householder Reflection)</span> $H^T H = I$ and $\\det(H) = -1$: preserves lengths and angles while flipping orientation across the dashed mirror line ($\\phi = ${param}^\\circ$). Also symmetric: $H^{-1} = H^T = H$!`;
        } else {
          verdictEl.innerHTML = `<span class="m10-pill bad">NOT Orthogonal (Despite $\\det(S) = +1$!)</span> Even though area is preserved ($\\det S = 1$), the unit circle stretches into an <strong>ellipse</strong> and lengths/angles distort whenever $k \\ne 0$. Exam reminder: $\\det(A) = \\pm 1$ does <em>not</em> guarantee orthogonality!`;
        }
        renderMath(verdictEl);

        drawIsoCanvas(M, Au, Av, rad);
      }

      function drawArrow(ctx, x0, y0, x1, y1, color, width, dashed) {
        ctx.save();
        ctx.strokeStyle = color;
        ctx.fillStyle = color;
        ctx.lineWidth = width;
        if (dashed) ctx.setLineDash([5, 4]);
        ctx.beginPath();
        ctx.moveTo(x0, y0);
        ctx.lineTo(x1, y1);
        ctx.stroke();
        ctx.setLineDash([]);
        const ang = Math.atan2(y1 - y0, x1 - x0);
        const hl = 9;
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x1 - hl * Math.cos(ang - 0.42), y1 - hl * Math.sin(ang - 0.42));
        ctx.lineTo(x1 - hl * Math.cos(ang + 0.42), y1 - hl * Math.sin(ang + 0.42));
        ctx.closePath();
        ctx.fill();
        ctx.restore();
      }

      function drawIsoCanvas(M, Au, Av, rad) {
        if (!isoCanvas.isConnected) return;
        const dpr = window.devicePixelRatio || 1;
        const rect = isoCanvas.getBoundingClientRect();
        const w = rect.width || 400, h = rect.height || 300;
        isoCanvas.width = Math.round(w * dpr);
        isoCanvas.height = Math.round(h * dpr);
        isoCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
        isoCtx.clearRect(0, 0, w, h);

        const { toPx, cx, cy, scale } = getIsoTransformCoords(w, h);

        // Grid
        isoCtx.strokeStyle = '#F1F5F9';
        isoCtx.lineWidth = 1;
        for (let g = -3; g <= 3; g++) {
          const [gx] = toPx([g, 0]);
          const [, gy] = toPx([0, g]);
          isoCtx.beginPath(); isoCtx.moveTo(gx, 0); isoCtx.lineTo(gx, h); isoCtx.stroke();
          isoCtx.beginPath(); isoCtx.moveTo(0, gy); isoCtx.lineTo(w, gy); isoCtx.stroke();
        }

        // Axes
        isoCtx.strokeStyle = '#CBD5E1';
        isoCtx.lineWidth = 1.4;
        isoCtx.beginPath(); isoCtx.moveTo(0, cy); isoCtx.lineTo(w, cy); isoCtx.stroke();
        isoCtx.beginPath(); isoCtx.moveTo(cx, 0); isoCtx.lineTo(cx, h); isoCtx.stroke();

        // Unit circle (dashed gray)
        isoCtx.save();
        isoCtx.strokeStyle = '#94A3B8';
        isoCtx.lineWidth = 1.4;
        isoCtx.setLineDash([4, 4]);
        isoCtx.beginPath();
        isoCtx.arc(cx, cy, scale, 0, Math.PI * 2);
        isoCtx.stroke();
        isoCtx.restore();

        // Image of unit circle under A (solid purple)
        isoCtx.beginPath();
        for (let i = 0; i <= 96; i++) {
          const t = (i / 96) * Math.PI * 2;
          const pt = matVec2(M, [Math.cos(t), Math.sin(t)]);
          const [px, py] = toPx(pt);
          if (i === 0) isoCtx.moveTo(px, py); else isoCtx.lineTo(px, py);
        }
        isoCtx.strokeStyle = 'rgba(124, 58, 237, 0.72)';
        isoCtx.lineWidth = 2.2;
        isoCtx.stroke();

        // Reflection mirror axis if in ref mode
        if (isoMode === 'ref') {
          const dx = Math.cos(rad) * 4, dy = Math.sin(rad) * 4;
          const [x1, y1] = toPx([-dx, -dy]);
          const [x2, y2] = toPx([dx, dy]);
          isoCtx.save();
          isoCtx.strokeStyle = '#DB2777';
          isoCtx.lineWidth = 1.8;
          isoCtx.setLineDash([6, 4]);
          isoCtx.beginPath(); isoCtx.moveTo(x1, y1); isoCtx.lineTo(x2, y2); isoCtx.stroke();
          isoCtx.restore();
        }

        // Draw original u, v (dashed lighter) and transformed Au, Av (solid bold)
        const [ux, uy] = toPx(uVec), [vx, vy] = toPx(vVec);
        const [Aux, Auy] = toPx(Au), [Avx, Avy] = toPx(Av);

        drawArrow(isoCtx, cx, cy, ux, uy, 'rgba(15,118,110,0.45)', 2, true);
        drawArrow(isoCtx, cx, cy, vx, vy, 'rgba(180,83,9,0.45)', 2, true);
        drawArrow(isoCtx, cx, cy, Aux, Auy, '#0F766E', 2.8, false);
        drawArrow(isoCtx, cx, cy, Avx, Avy, '#B45309', 2.8, false);

        // Draggable handles on u and v
        [{ x: ux, y: uy, col: '#0F766E', lbl: 'u' }, { x: vx, y: vy, col: '#B45309', lbl: 'v' }].forEach(p => {
          isoCtx.beginPath();
          isoCtx.arc(p.x, p.y, 6, 0, Math.PI * 2);
          isoCtx.fillStyle = '#FFFFFF';
          isoCtx.fill();
          isoCtx.strokeStyle = p.col;
          isoCtx.lineWidth = 2.2;
          isoCtx.stroke();
          isoCtx.fillStyle = p.col;
          isoCtx.font = 'bold 12px sans-serif';
          isoCtx.fillText(p.lbl, p.x + 8, p.y - 6);
        });

        isoCtx.font = 'bold 12px sans-serif';
        isoCtx.fillStyle = '#0F766E'; isoCtx.fillText('Au', Aux + 7, Auy + 4);
        isoCtx.fillStyle = '#B45309'; isoCtx.fillText('Av', Avx + 7, Avy + 4);
      }

      // Pointer dragging for u and v on Lab 1 canvas
      function getCanvasPoint(evt) {
        const rect = isoCanvas.getBoundingClientRect();
        const clientX = evt.touches ? evt.touches[0].clientX : evt.clientX;
        const clientY = evt.touches ? evt.touches[0].clientY : evt.clientY;
        return { x: clientX - rect.left, y: clientY - rect.top, w: rect.width, h: rect.height };
      }

      function onIsoDown(evt) {
        const { x, y, w, h } = getCanvasPoint(evt);
        const { toPx } = getIsoTransformCoords(w, h);
        const [ux, uy] = toPx(uVec), [vx, vy] = toPx(vVec);
        const du = Math.hypot(x - ux, y - uy), dv = Math.hypot(x - vx, y - vy);
        if (du < 18 || dv < 18) {
          dragging = du <= dv ? 'u' : 'v';
          isoCanvas.style.cursor = 'grabbing';
          evt.preventDefault();
        }
      }

      function onIsoMove(evt) {
        if (!dragging) return;
        const { x, y, w, h } = getCanvasPoint(evt);
        const { toMath } = getIsoTransformCoords(w, h);
        const [mx, my] = toMath(x, y);
        const clamped = [
          Math.max(-3.2, Math.min(3.2, Math.round(mx * 10) / 10)),
          Math.max(-3.2, Math.min(3.2, Math.round(my * 10) / 10))
        ];
        if (dragging === 'u') uVec = clamped; else vVec = clamped;
        updateIsoLab();
        evt.preventDefault();
      }

      function onIsoUp() {
        dragging = null;
        isoCanvas.style.cursor = 'grab';
      }

      isoCanvas.addEventListener('mousedown', onIsoDown);
      window.addEventListener('mousemove', onIsoMove);
      window.addEventListener('mouseup', onIsoUp);
      isoCanvas.addEventListener('touchstart', onIsoDown, { passive: false });
      window.addEventListener('touchmove', onIsoMove, { passive: false });
      window.addEventListener('touchend', onIsoUp);

      isoSlider.addEventListener('input', updateIsoLab);
      container.querySelectorAll('[data-iso]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-iso]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          isoMode = btn.dataset.iso;
          if (isoMode === 'rot') {
            isoSlider.min = -180; isoSlider.max = 180; isoSlider.step = 1; isoSlider.value = 50;
            isoLabel.innerHTML = 'Rotation Angle $\\theta$:';
          } else if (isoMode === 'ref') {
            isoSlider.min = -90; isoSlider.max = 90; isoSlider.step = 1; isoSlider.value = 30;
            isoLabel.innerHTML = 'Mirror Axis Angle $\\phi$:';
          } else {
            isoSlider.min = -1.5; isoSlider.max = 1.5; isoSlider.step = 0.05; isoSlider.value = 0.8;
            isoLabel.innerHTML = 'Shear Factor $k$:';
          }
          renderMath(isoLabel);
          updateIsoLab();
        });
      });

      // =====================================================================
      // LAB 2: L^2[-1,1] & LEGENDRE ORTHOGONAL PROJECTION
      // =====================================================================
      let l2Key = 'runge';
      const degSlider = container.querySelector('#m10DegSlider');
      const degVal = container.querySelector('#m10DegVal');
      const showTaylorChk = container.querySelector('#m10ShowTaylor');
      const l2Canvas = container.querySelector('#m10L2Canvas');
      const l2Ctx = l2Canvas.getContext('2d');

      function updateL2Lab() {
        const n = parseInt(degSlider.value, 10);
        degVal.textContent = `n = ${n}`;
        const spec = L2_FUNCS[l2Key];
        const f = spec.fn;

        // Compute Legendre coefficients a_k = <f, P_k> / (2/(2k+1))
        const a = [];
        for (let k = 0; k <= n; k++) {
          const ip = innerL2(f, x => legendreP(k, x));
          a.push(ip * ((2 * k + 1) / 2));
        }

        const fHat = x => {
          let s = 0;
          for (let k = 0; k <= n; k++) s += a[k] * legendreP(k, x);
          return s;
        };

        const fNormSq = innerL2(f, f);
        let projNormSq = 0;
        for (let k = 0; k <= n; k++) projNormSq += (a[k] * a[k] * 2) / (2 * k + 1);
        const l2Err = Math.sqrt(Math.max(0, innerL2(x => f(x) - fHat(x), x => f(x) - fHat(x))));
        const energyPct = fNormSq > 0 ? Math.min(100, (projNormSq / fNormSq) * 100) : 100;

        let taylorErrStr = 'N/A';
        if (spec.hasTaylor) {
          const tErr = Math.sqrt(Math.max(0, innerL2(x => f(x) - spec.taylor(n, x), x => f(x) - spec.taylor(n, x))));
          taylorErrStr = tErr > 99 ? tErr.toExponential(1) : clean(tErr, 4);
        }

        const statsEl = container.querySelector('#m10L2Stats');
        statsEl.innerHTML = `
          <div class="m10-stat-card">
            <div class="m10-stat-lbl">Legendre $L^2$ Error</div>
            <div class="m10-stat-val" style="color:#0F766E;">${clean(l2Err, 4)}</div>
          </div>
          <div class="m10-stat-card">
            <div class="m10-stat-lbl">Taylor $L^2$ Error</div>
            <div class="m10-stat-val" style="color:${spec.hasTaylor ? '#DC2626' : '#64748B'};">${taylorErrStr}</div>
          </div>
          <div class="m10-stat-card">
            <div class="m10-stat-lbl">Parseval Energy Captured</div>
            <div class="m10-stat-val" style="color:#1D4ED8;">${energyPct.toFixed(2)}%</div>
          </div>
          <div class="m10-stat-card">
            <div class="m10-stat-lbl">$\\|f\\|_{L^2}^2$ vs $\\|\\hat{f}_n\\|_{L^2}^2$</div>
            <div class="m10-stat-val" style="font-size:0.84rem;">${clean(fNormSq, 3)} / ${clean(projNormSq, 3)}</div>
          </div>
        `;
        renderMath(statsEl);

        // Format compact expansion terms
        const terms = a
          .map((coef, k) => (Math.abs(coef) > 1e-4 ? `${coef >= 0 && k > 0 ? '+' : ''}${clean(coef, 3)}P_${k}(x)` : ''))
          .filter(Boolean)
          .join(' ');
        const readoutEl = container.querySelector('#m10L2Readout');
        readoutEl.innerHTML = `
          <div><strong>Orthogonal Projection onto $\\mathcal{P}_${n}$:</strong> $\\hat{f}_${n}(x) \\approx ${terms || '0'}$</div>
          <div style="margin-top:5px; color:#334155;">${spec.note}</div>
        `;
        renderMath(readoutEl);

        drawL2Canvas(spec, fHat, n);
      }

      function drawL2Canvas(spec, fHat, n) {
        if (!l2Canvas.isConnected) return;
        const dpr = window.devicePixelRatio || 1;
        const rect = l2Canvas.getBoundingClientRect();
        const w = rect.width || 680, h = rect.height || 280;
        l2Canvas.width = Math.round(w * dpr);
        l2Canvas.height = Math.round(h * dpr);
        l2Ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        l2Ctx.clearRect(0, 0, w, h);

        const padL = 42, padR = 16, padT = 26, padB = 26;
        const plotW = w - padL - padR, plotH = h - padT - padB;
        const SAMPLES = 240;

        // Compute y bounds from f and fHat (so exploding Taylor doesn't squash the plot)
        let yMin = Infinity, yMax = -Infinity;
        for (let i = 0; i <= SAMPLES; i++) {
          const x = -1 + (2 * i) / SAMPLES;
          const yf = spec.fn(x), yh = fHat(x);
          yMin = Math.min(yMin, yf, yh);
          yMax = Math.max(yMax, yf, yh);
        }
        const padY = Math.max(0.25, (yMax - yMin) * 0.16);
        yMin -= padY; yMax += padY;

        const mapX = x => padL + ((x + 1) / 2) * plotW;
        const mapY = y => padT + plotH - ((y - yMin) / (yMax - yMin || 1)) * plotH;

        // Axes
        l2Ctx.strokeStyle = '#E2E8F0';
        l2Ctx.lineWidth = 1.2;
        if (yMin < 0 && yMax > 0) {
          l2Ctx.beginPath(); l2Ctx.moveTo(padL, mapY(0)); l2Ctx.lineTo(padL + plotW, mapY(0)); l2Ctx.stroke();
        }
        l2Ctx.beginPath(); l2Ctx.moveTo(mapX(0), padT); l2Ctx.lineTo(mapX(0), padT + plotH); l2Ctx.stroke();

        // X ticks
        l2Ctx.fillStyle = '#64748B';
        l2Ctx.font = '10.5px monospace';
        l2Ctx.textAlign = 'center';
        [-1, -0.5, 0, 0.5, 1].forEach(xv => {
          l2Ctx.fillText(String(xv), mapX(xv), padT + plotH + 16);
        });

        // Clip to plot area for curves
        l2Ctx.save();
        l2Ctx.beginPath();
        l2Ctx.rect(padL, padT, plotW, plotH);
        l2Ctx.clip();

        // 1) Optional Taylor polynomial (dashed red)
        if (showTaylorChk.checked && spec.hasTaylor) {
          l2Ctx.beginPath();
          for (let i = 0; i <= SAMPLES; i++) {
            const x = -1 + (2 * i) / SAMPLES;
            const yt = Math.max(yMin - 5, Math.min(yMax + 5, spec.taylor(n, x)));
            if (i === 0) l2Ctx.moveTo(mapX(x), mapY(yt)); else l2Ctx.lineTo(mapX(x), mapY(yt));
          }
          l2Ctx.strokeStyle = 'rgba(220, 38, 38, 0.8)';
          l2Ctx.lineWidth = 2;
          l2Ctx.setLineDash([6, 4]);
          l2Ctx.stroke();
          l2Ctx.setLineDash([]);
        }

        // 2) True function f(x) (dark slate)
        l2Ctx.beginPath();
        for (let i = 0; i <= SAMPLES; i++) {
          const x = -1 + (2 * i) / SAMPLES;
          const yf = spec.fn(x);
          if (i === 0) l2Ctx.moveTo(mapX(x), mapY(yf)); else l2Ctx.lineTo(mapX(x), mapY(yf));
        }
        l2Ctx.strokeStyle = '#334155';
        l2Ctx.lineWidth = 2.2;
        l2Ctx.stroke();

        // 3) Legendre L^2 Orthogonal Projection f_hat_n(x) (bold teal)
        l2Ctx.beginPath();
        for (let i = 0; i <= SAMPLES; i++) {
          const x = -1 + (2 * i) / SAMPLES;
          const yh = fHat(x);
          if (i === 0) l2Ctx.moveTo(mapX(x), mapY(yh)); else l2Ctx.lineTo(mapX(x), mapY(yh));
        }
        l2Ctx.strokeStyle = '#0F766E';
        l2Ctx.lineWidth = 2.8;
        l2Ctx.stroke();
        l2Ctx.restore();

        // Legend
        l2Ctx.font = 'bold 11px sans-serif';
        l2Ctx.textAlign = 'left';
        l2Ctx.fillStyle = '#334155'; l2Ctx.fillText('━ Target f(x)', padL, 16);
        l2Ctx.fillStyle = '#0F766E'; l2Ctx.fillText(`━ L² Legendre Projection (deg ${n})`, padL + 100, 16);
        if (showTaylorChk.checked && spec.hasTaylor) {
          l2Ctx.fillStyle = '#DC2626'; l2Ctx.fillText(`┅ Taylor T_${n}(x)`, padL + 305, 16);
        }
      }

      degSlider.addEventListener('input', updateL2Lab);
      showTaylorChk.addEventListener('change', updateL2Lab);
      container.querySelectorAll('[data-l2]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-l2]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          l2Key = btn.dataset.l2;
          updateL2Lab();
        });
      });

      // =====================================================================
      // SECTION C: TOPIC 6 COMPREHENSIVE EXAM SIMULATOR
      // =====================================================================
      let examOrder = EXAM_BANK.slice();
      let examMode = 'practice'; // 'practice' | 'exam'
      let curIdx = 0;
      let userAnswers = new Array(EXAM_BANK.length).fill(null);
      let examGraded = false;

      // Timer state
      let timerSeconds = 0;
      let timerId = null;
      const timerDisplay = container.querySelector('#m10TimerDisplay');
      const timerToggleBtn = container.querySelector('#m10TimerToggleBtn');

      function fmtTime(s) {
        const m = String(Math.floor(s / 60)).padStart(2, '0');
        const sec = String(s % 60).padStart(2, '0');
        return `${m}:${sec}`;
      }

      function stopTimer() {
        if (timerId) { clearInterval(timerId); timerId = null; }
        timerToggleBtn.textContent = '▶ Start Timer';
      }

      timerToggleBtn.addEventListener('click', () => {
        if (timerId) {
          stopTimer();
        } else {
          timerToggleBtn.textContent = '⏸ Pause';
          timerId = setInterval(() => {
            if (!timerDisplay.isConnected) { clearInterval(timerId); timerId = null; return; }
            timerSeconds++;
            timerDisplay.textContent = fmtTime(timerSeconds);
          }, 1000);
        }
      });

      function renderExamUI() {
        // Question navigation bar
        const navEl = container.querySelector('#m10QNav');
        navEl.innerHTML = examOrder.map((qObj, idx) => {
          const picked = userAnswers[idx];
          let cls = 'm10-qdot';
          if (idx === curIdx) cls += ' active';
          const showResult = examGraded || (examMode === 'practice' && picked !== null);
          if (showResult && picked !== null) {
            cls += picked === qObj.ans ? ' correct' : ' wrong';
          } else if (picked !== null) {
            cls += ' answered';
          }
          return `<button class="${cls}" data-qjump="${idx}">Q${idx + 1}</button>`;
        }).join('');

        navEl.querySelectorAll('[data-qjump]').forEach(btn => {
          btn.addEventListener('click', () => {
            curIdx = parseInt(btn.dataset.qjump, 10);
            renderExamUI();
          });
        });

        // Active question panel
        const qObj = examOrder[curIdx];
        const picked = userAnswers[curIdx];
        const showFeedback = examGraded || (examMode === 'practice' && picked !== null);
        const qBox = container.querySelector('#m10QuestionBox');

        qBox.innerHTML = `
          <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:8px;">
            <span class="m10-pill info">${qObj.pillar}</span>
            <span style="font-size:0.8rem; color:var(--text-muted); font-weight:700;">Question ${curIdx + 1} of ${examOrder.length}</span>
          </div>
          <div style="font-size:0.95rem; font-weight:600; color:#0F172A; margin-bottom:10px;">
            ${qObj.q}
          </div>
          <div>
            ${qObj.opts.map((opt, oi) => {
              let cls = 'm10-opt-btn';
              if (showFeedback) {
                if (oi === qObj.ans) cls += ' correct';
                else if (oi === picked) cls += ' wrong';
              } else if (oi === picked) {
                cls += ' selected';
              }
              return `<button class="${cls}" data-eopt="${oi}"><strong>${String.fromCharCode(65 + oi)}.</strong> ${opt}</button>`;
            }).join('')}
          </div>
          ${showFeedback ? `
            <div style="margin-top:12px; padding:10px 12px; border-radius:8px; background:${picked === qObj.ans ? '#F0FDF4' : '#FEF2F2'}; border:1px solid ${picked === qObj.ans ? '#BBF7D0' : '#FECACA'}; font-size:0.86rem;">
              <strong style="color:${picked === qObj.ans ? '#166534' : '#991B1B'};">
                ${picked === null ? 'Unanswered' : (picked === qObj.ans ? '✓ Correct!' : '✗ Incorrect')}
              </strong> — ${qObj.exp}
            </div>
          ` : ''}
        `;

        qBox.querySelectorAll('[data-eopt]').forEach(btn => {
          btn.addEventListener('click', () => {
            if (examGraded) return;
            userAnswers[curIdx] = parseInt(btn.dataset.eopt, 10);
            renderExamUI();
          });
        });

        renderMath(qBox);
      }

      function renderMasteryReport() {
        examGraded = true;
        stopTimer();
        renderExamUI();

        let totalCorrect = 0;
        const pillarStats = {};
        examOrder.forEach((qObj, idx) => {
          if (!pillarStats[qObj.pillar]) pillarStats[qObj.pillar] = { correct: 0, total: 0 };
          pillarStats[qObj.pillar].total++;
          if (userAnswers[idx] === qObj.ans) {
            totalCorrect++;
            pillarStats[qObj.pillar].correct++;
          }
        });

        const pct = Math.round((totalCorrect / examOrder.length) * 100);
        const tier = pct >= 85 ? 'Exam Ready — Distinction Level!' : (pct >= 65 ? 'Solid Foundation — Polish Weak Pillars' : 'Review Recommended Before Exam');
        const tierCol = pct >= 85 ? '#166534' : (pct >= 65 ? '#B45309' : '#991B1B');

        const reportEl = container.querySelector('#m10ReportBox');
        reportEl.style.display = 'block';
        reportEl.innerHTML = `
          <div class="m10-panel" style="border-color:#CBD5E1; background:#F8FAFC;">
            <div style="display:flex; flex-wrap:wrap; justify-content:space-between; align-items:center; gap:8px; margin-bottom:10px;">
              <div>
                <div style="font-size:1.05rem; font-weight:800; color:${tierCol};">Score: ${totalCorrect} / ${examOrder.length} (${pct}%) — ${tier}</div>
                <div style="font-size:0.82rem; color:var(--text-muted);">Elapsed time: ${fmtTime(timerSeconds)} · Click any question pill Q1–Q10 above to review full solutions.</div>
              </div>
            </div>
            <div style="display:grid; grid-template-columns:repeat(auto-fit, minmax(220px, 1fr)); gap:10px;">
              ${Object.entries(pillarStats).map(([pillar, st]) => {
                const pPct = Math.round((st.correct / st.total) * 100);
                const barCol = pPct === 100 ? '#16A34A' : (pPct >= 50 ? '#D97706' : '#DC2626');
                return `
                  <div style="background:#FFFFFF; border:1px solid var(--border); border-radius:8px; padding:9px 11px;">
                    <div style="display:flex; justify-content:space-between; font-size:0.8rem; font-weight:700;">
                      <span>${pillar}</span>
                      <span style="color:${barCol};">${st.correct}/${st.total}</span>
                    </div>
                    <div class="m10-bar-track">
                      <div class="m10-bar-fill" style="width:${pPct}%; background:${barCol};"></div>
                    </div>
                  </div>
                `;
              }).join('')}
            </div>
          </div>
        `;
      }

      function resetExam(shuffle) {
        if (shuffle) {
          examOrder = EXAM_BANK.slice().sort(() => Math.random() - 0.5);
        }
        curIdx = 0;
        userAnswers = new Array(examOrder.length).fill(null);
        examGraded = false;
        stopTimer();
        timerSeconds = 0;
        timerDisplay.textContent = '00:00';
        container.querySelector('#m10ReportBox').style.display = 'none';
        renderExamUI();
      }

      container.querySelectorAll('[data-emode]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-emode]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          examMode = btn.dataset.emode;
          resetExam(false);
        });
      });

      container.querySelector('#m10ShuffleBtn').addEventListener('click', () => resetExam(true));
      container.querySelector('#m10ResetExamBtn').addEventListener('click', () => resetExam(false));
      container.querySelector('#m10PrevQ').addEventListener('click', () => {
        curIdx = (curIdx - 1 + examOrder.length) % examOrder.length;
        renderExamUI();
      });
      container.querySelector('#m10NextQ').addEventListener('click', () => {
        curIdx = (curIdx + 1) % examOrder.length;
        renderExamUI();
      });
      container.querySelector('#m10SubmitExamBtn').addEventListener('click', renderMasteryReport);

      // Initial renders
      updateIsoLab();
      updateL2Lab();
      renderExamUI();

      window.addEventListener('resize', () => {
        updateIsoLab();
        updateL2Lab();
      });
    }
  });
})();
