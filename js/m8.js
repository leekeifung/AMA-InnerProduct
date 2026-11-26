/**
 * AMA1751 Linear Algebra — Topic 6: Inner Product
 * Module 8: Comprehensive Mastery Arena, Procedural Problem Engine & Mock Exam Simulator
 */
(function () {
  function clean(x, d = 4) {
    if (!Number.isFinite(x)) return 'NaN';
    if (Math.abs(x) < 1e-10) return '0';
    return String(Number(x.toFixed(d)));
  }

  function pick(arr) {
    return arr[Math.floor(Math.random() * arr.length)];
  }

  function shuffle(arr) {
    const a = arr.slice();
    for (let i = a.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [a[i], a[j]] = [a[j], a[i]];
    }
    return a;
  }

  /** Parse integers, decimals, or exact fractions like "-3/4" */
  function parseNum(str) {
    const s = String(str || '').trim().replace(/\s+/g, '');
    if (!s) return NaN;
    if (/^[+-]?\d+(\.\d+)?\/[+-]?\d+(\.\d+)?$/.test(s)) {
      const [num, den] = s.split('/').map(Number);
      return den === 0 ? NaN : num / den;
    }
    return Number(s);
  }

  function colVecTex(v) {
    return `\\begin{bmatrix} ${v.map(x => clean(x, 3)).join(' \\\\ ')} \\end{bmatrix}`;
  }

  function matTex(M) {
    return `\\begin{bmatrix} ${M.map(r => r.map(x => clean(x, 3)).join(' & ')).join(' \\\\ ')} \\end{bmatrix}`;
  }

  // =========================================================================
  // FORMULA ATLAS DATA
  // =========================================================================
  const FORMULAS = [
    {
      cat: 'norm',
      tag: 'Norm & Orthogonality',
      title: 'Dot Product, Norm, Angle & Cauchy–Schwarz',
      body: `$$\\|\\mathbf{u}\\| = \\sqrt{\\mathbf{u}\\cdot\\mathbf{u}}, \\quad \\cos\\theta = \\frac{\\mathbf{u}\\cdot\\mathbf{v}}{\\|\\mathbf{u}\\|\\|\\mathbf{v}\\|}, \\quad |\\mathbf{u}\\cdot\\mathbf{v}| \\le \\|\\mathbf{u}\\|\\|\\mathbf{v}\\|$$
             <strong>Pythagoras:</strong> $\\|\\mathbf{u}+\\mathbf{v}\\|^2 = \\|\\mathbf{u}\\|^2 + \\|\\mathbf{v}\\|^2 \\iff \\mathbf{u} \\perp \\mathbf{v}$.`
    },
    {
      cat: 'proj',
      tag: 'Projections & Subspaces',
      title: 'Orthogonal Complements & 4 Fundamental Subspaces',
      body: `$$(\\operatorname{Row} A)^\\perp = \\operatorname{Nul} A, \\qquad (\\operatorname{Col} A)^\\perp = \\operatorname{Nul} A^T$$
             <strong>Dimension Law:</strong> $\\dim W + \\dim W^\\perp = n$ in $\\mathbb{R}^n$. Every $\\mathbf{y}\\in\\mathbb{R}^n$ splits uniquely as $\\mathbf{y} = \\hat{\\mathbf{y}} + \\mathbf{z}$ ($\\hat{\\mathbf{y}}\\in W,\\; \\mathbf{z}\\in W^\\perp$).`
    },
    {
      cat: 'proj',
      tag: 'Projections & Subspaces',
      title: 'Orthogonal Projection & Best Approximation',
      body: `For an orthogonal basis $\\{\\mathbf{u}_1,\\dots,\\mathbf{u}_k\\}$ of $W$:
             $$\\operatorname{proj}_W(\\mathbf{y}) = \\sum_{i=1}^k \\frac{\\mathbf{y}\\cdot\\mathbf{u}_i}{\\mathbf{u}_i\\cdot\\mathbf{u}_i}\\mathbf{u}_i = U U^T \\mathbf{y} \\quad (U \\text{ orthonormal columns})$$
             Minimizes $\\|\\mathbf{y} - \\mathbf{w}\\|$ over all $\\mathbf{w} \\in W$; distance $= \\|\\mathbf{y} - \\hat{\\mathbf{y}}\\|$.`
    },
    {
      cat: 'gsqr',
      tag: 'Gram–Schmidt & QR',
      title: 'Gram–Schmidt Process & Reduced QR Factorization',
      body: `$$\\mathbf{v}_k = \\mathbf{x}_k - \\sum_{j=1}^{k-1} \\frac{\\mathbf{x}_k\\cdot\\mathbf{v}_j}{\\mathbf{v}_j\\cdot\\mathbf{v}_j}\\mathbf{v}_j, \\qquad \\mathbf{q}_k = \\frac{\\mathbf{v}_k}{\\|\\mathbf{v}_k\\|}$$
             Produces $A = QR$ ($Q^T Q = I_n$, $R = Q^T A$ upper triangular with $r_{ii} > 0$).`
    },
    {
      cat: 'lsq',
      tag: 'Least Squares',
      title: 'Normal Equations & QR Least-Squares Solver',
      body: `$$A^T A \\hat{\\mathbf{x}} = A^T \\mathbf{b} \\quad \\Longleftrightarrow \\quad R\\hat{\\mathbf{x}} = Q^T \\mathbf{b}$$
             Unique solution $\\iff$ columns of $A$ are linearly independent $\\iff$ $A^T A$ is invertible. Residual $\\mathbf{e} = \\mathbf{b} - A\\hat{\\mathbf{x}} \\perp \\operatorname{Col} A$.`
    },
    {
      cat: 'ips',
      tag: 'Inner Product Spaces & Fourier',
      title: 'General Inner Product, Legendre & Fourier Series',
      body: `<strong>4 Axioms:</strong> Symmetry, Additivity, Homogeneity, Positivity ($\\langle \\mathbf{u},\\mathbf{u}\\rangle \\ge 0$, $=0\\iff \\mathbf{u}=\\mathbf{0}$).<br/>
             <strong>$L^2$ Projection:</strong> $\\hat{f} = \\sum \\frac{\\langle f, \\phi_k\\rangle}{\\langle \\phi_k, \\phi_k\\rangle}\\phi_k$ (Legendre on $[-1,1]$, Fourier $\\{1,\\cos kt, \\sin kt\\}$ on $[-\\pi,\\pi]$).`
    }
  ];

  // =========================================================================
  // PROCEDURAL PROBLEM ENGINE (EXACT INTEGER SOLUTIONS)
  // =========================================================================
  function genProjectionProblem() {
    const pairs = [
      { u1: [1, 1, 0],  u2: [1, -1, 1], n: [1, -1, -2] },
      { u1: [1, 0, 1],  u2: [1, 2, -1], n: [-1, 1, 1] },
      { u1: [2, 1, 0],  u2: [-1, 2, 0], n: [0, 0, 1] },
      { u1: [1, -1, 0], u2: [1, 1, 2],  n: [-1, -1, 1] }
    ];
    const p = pick(pairs);
    const c1 = pick([-2, -1, 1, 2]);
    const c2 = pick([-2, -1, 1, 2]);
    const c3 = pick([-2, -1, 1, 2]);
    const yhat = [0, 1, 2].map(i => c1 * p.u1[i] + c2 * p.u2[i]);
    const y = [0, 1, 2].map(i => yhat[i] + c3 * p.n[i]);
    const d1 = p.u1.reduce((s, x) => s + x * x, 0);
    const d2 = p.u2.reduce((s, x) => s + x * x, 0);
    const yu1 = y.reduce((s, x, i) => s + x * p.u1[i], 0);
    const yu2 = y.reduce((s, x, i) => s + x * p.u2[i], 0);
    const errSq = [0, 1, 2].reduce((s, i) => s + Math.pow(y[i] - yhat[i], 2), 0);

    return {
      title: 'Orthogonal Projection onto a 2D Subspace',
      prompt: `Let $W = \\operatorname{span}\\{\\mathbf{u}_1, \\mathbf{u}_2\\}$ with orthogonal basis $\\mathbf{u}_1 = ${colVecTex(p.u1)}$, $\\mathbf{u}_2 = ${colVecTex(p.u2)}$. Find the orthogonal projection $\\hat{\\mathbf{y}} = \\operatorname{proj}_W(\\mathbf{y})$ of $\\mathbf{y} = ${colVecTex(y)}$ onto $W$.`,
      fields: [
        { id: 'y1', label: 'ŷ₁', ans: yhat[0] },
        { id: 'y2', label: 'ŷ₂', ans: yhat[1] },
        { id: 'y3', label: 'ŷ₃', ans: yhat[2] }
      ],
      solution: `
        <strong>Step 1 (Verify orthogonality):</strong> $\\mathbf{u}_1\\cdot\\mathbf{u}_2 = 0$, $\\|\\mathbf{u}_1\\|^2 = ${d1}$, $\\|\\mathbf{u}_2\\|^2 = ${d2}$.<br/>
        <strong>Step 2 (Compute coordinates):</strong>
        $$c_1 = \\frac{\\mathbf{y}\\cdot\\mathbf{u}_1}{\\mathbf{u}_1\\cdot\\mathbf{u}_1} = \\frac{${yu1}}{${d1}} = ${c1}, \\qquad c_2 = \\frac{\\mathbf{y}\\cdot\\mathbf{u}_2}{\\mathbf{u}_2\\cdot\\mathbf{u}_2} = \\frac{${yu2}}{${d2}} = ${c2}$$
        <strong>Step 3 (Assemble projection):</strong>
        $$\\hat{\\mathbf{y}} = (${c1})${colVecTex(p.u1)} + (${c2})${colVecTex(p.u2)} = ${colVecTex(yhat)}$$
        <em>Bonus check:</em> Residual $\\mathbf{y} - \\hat{\\mathbf{y}} = ${colVecTex(y.map((v, i) => v - yhat[i]))}$ is orthogonal to both $\\mathbf{u}_1, \\mathbf{u}_2$, with squared distance $\\|\\mathbf{y} - \\hat{\\mathbf{y}}\\|^2 = ${errSq}$.
      `
    };
  }

  function genGramSchmidtProblem() {
    const pairs = [
      { u1: [1, 1, 0],  u2: [1, -1, 1] },
      { u1: [1, 0, 1],  u2: [1, 2, -1] },
      { u1: [1, -1, 1], u2: [1, 1, 0] },
      { u1: [2, 1, 0],  u2: [-1, 2, 1] },
      { u1: [0, 1, 2],  u2: [1, -2, 1] }
    ];
    const p = pick(pairs);
    const k = pick([-2, -1, 1, 2]);
    const s = pick([-1, 1]);
    const v1 = p.u1.slice();
    const v2 = p.u2.map(x => s * x);
    const x2 = [0, 1, 2].map(i => k * v1[i] + v2[i]);
    const d1 = v1.reduce((acc, x) => acc + x * x, 0);
    const x2v1 = x2.reduce((acc, x, i) => acc + x * v1[i], 0);

    return {
      title: 'Gram–Schmidt Orthogonalization Step',
      prompt: `Apply the Gram–Schmidt process to $\\mathbf{x}_1 = ${colVecTex(v1)}$ and $\\mathbf{x}_2 = ${colVecTex(x2)}$ (setting $\\mathbf{v}_1 = \\mathbf{x}_1$). Compute the exact orthogonal vector $\\mathbf{v}_2 = \\mathbf{x}_2 - \\operatorname{proj}_{\\mathbf{v}_1}(\\mathbf{x}_2)$.`,
      fields: [
        { id: 'v1', label: 'v₂,₁', ans: v2[0] },
        { id: 'v2', label: 'v₂,₂', ans: v2[1] },
        { id: 'v3', label: 'v₂,₃', ans: v2[2] }
      ],
      solution: `
        <strong>Step 1 (Projection scalar):</strong> $\\mathbf{x}_2\\cdot\\mathbf{v}_1 = ${x2v1}$ and $\\mathbf{v}_1\\cdot\\mathbf{v}_1 = ${d1}$, so
        $$\\frac{\\mathbf{x}_2\\cdot\\mathbf{v}_1}{\\mathbf{v}_1\\cdot\\mathbf{v}_1} = \\frac{${x2v1}}{${d1}} = ${k}$$
        <strong>Step 2 (Subtract projection):</strong>
        $$\\mathbf{v}_2 = ${colVecTex(x2)} - (${k})${colVecTex(v1)} = ${colVecTex(v2)}$$
        <strong>Check:</strong> $\\mathbf{v}_1\\cdot\\mathbf{v}_2 = 0$ ✓.
      `
    };
  }

  function genLeastSquaresProblem() {
    const mats = [
      { A: [[1, 0], [1, 1], [0, 1]], n: [1, -1, 1] },
      { A: [[1, 1], [1, -1], [1, 0]], n: [1, 1, -2] },
      { A: [[1, 0], [0, 1], [1, 1]], n: [-1, -1, 1] },
      { A: [[2, 1], [1, 0], [0, 1]], n: [1, -2, -1] }
    ];
    const m = pick(mats);
    const x1 = pick([-2, -1, 1, 2, 3]);
    const x2 = pick([-2, -1, 1, 2, 3]);
    const c = pick([-2, -1, 1, 2]);
    const b = [0, 1, 2].map(i => m.A[i][0] * x1 + m.A[i][1] * x2 + c * m.n[i]);

    // Compute A^T A and A^T b
    const ata = [
      [
        [0, 1, 2].reduce((s, i) => s + m.A[i][0] * m.A[i][0], 0),
        [0, 1, 2].reduce((s, i) => s + m.A[i][0] * m.A[i][1], 0)
      ],
      [
        [0, 1, 2].reduce((s, i) => s + m.A[i][1] * m.A[i][0], 0),
        [0, 1, 2].reduce((s, i) => s + m.A[i][1] * m.A[i][1], 0)
      ]
    ];
    const atb = [
      [0, 1, 2].reduce((s, i) => s + m.A[i][0] * b[i], 0),
      [0, 1, 2].reduce((s, i) => s + m.A[i][1] * b[i], 0)
    ];

    return {
      title: 'Least-Squares Solution via Normal Equations',
      prompt: `Find the unique least-squares solution $\\hat{\\mathbf{x}} = \\begin{bmatrix} \\hat{x}_1 \\\\ \\hat{x}_2 \\end{bmatrix}$ to the inconsistent system $A\\mathbf{x} = \\mathbf{b}$, where $A = ${matTex(m.A)}$ and $\\mathbf{b} = ${colVecTex(b)}$.`,
      fields: [
        { id: 'x1', label: 'x̂₁', ans: x1 },
        { id: 'x2', label: 'x̂₂', ans: x2 }
      ],
      solution: `
        <strong>Step 1 (Form Normal Equations $A^T A \\hat{\\mathbf{x}} = A^T \\mathbf{b}$):</strong>
        $$A^T A = ${matTex(ata)}, \\qquad A^T \\mathbf{b} = ${colVecTex(atb)}$$
        <strong>Step 2 (Solve the $2\\times 2$ system):</strong>
        $$${matTex(ata)}\\begin{bmatrix}\\hat{x}_1\\\\\\hat{x}_2\\end{bmatrix} = ${colVecTex(atb)} \\implies \\hat{\\mathbf{x}} = ${colVecTex([x1, x2])}$$
        <strong>Step 3 (Verify residual orthogonality):</strong> $\\mathbf{e} = \\mathbf{b} - A\\hat{\\mathbf{x}} = ${colVecTex(m.n.map(v => c * v))}$ satisfies $A^T \\mathbf{e} = \\mathbf{0}$ ✓.
      `
    };
  }

  function genQRProblem() {
    const qSet = [
      { q1: [2, 2, 1], q2: [-1, 2, -2], q3: [2, -1, -2], den: 3 },
      { q1: [1, 2, 2], q2: [2, -2, 1],  q3: [2, 1, -2],  den: 3 }
    ];
    const q = pick(qSet);
    const x1 = pick([-2, -1, 1, 2]);
    const x2 = pick([-2, -1, 1, 2]);
    const s = pick([-1, 1]);
    const r12 = 3 * s;
    const R = [[3, r12], [0, 3]];
    const k1 = x1 + s * x2;
    const k2 = x2;
    const c = pick([-1, 1]);
    const b = [0, 1, 2].map(i => k1 * q.q1[i] + k2 * q.q2[i] + c * q.q3[i]);
    const qtb = [3 * k1, 3 * k2];

    const Qtex = `\\frac{1}{3}${matTex([[q.q1[0], q.q2[0]], [q.q1[1], q.q2[1]], [q.q1[2], q.q2[2]]])}`;

    return {
      title: 'Least-Squares via QR Factorization (R x̂ = Qᵀ b)',
      prompt: `Given the reduced factorization $A = QR$ with $Q = ${Qtex}$ and $R = ${matTex(R)}$, use $R\\hat{\\mathbf{x}} = Q^T \\mathbf{b}$ to find the least-squares solution $\\hat{\\mathbf{x}}$ for $\\mathbf{b} = ${colVecTex(b)}$.`,
      fields: [
        { id: 'x1', label: 'x̂₁', ans: x1 },
        { id: 'x2', label: 'x̂₂', ans: x2 }
      ],
      solution: `
        <strong>Step 1 (Compute $Q^T \\mathbf{b}$):</strong>
        $$Q^T \\mathbf{b} = \\frac{1}{3}\\begin{bmatrix} ${q.q1.join(' & ')} \\\\ ${q.q2.join(' & ')} \\end{bmatrix}${colVecTex(b)} = ${colVecTex(qtb)}$$
        <strong>Step 2 (Back-substitute in $R\\hat{\\mathbf{x}} = Q^T \\mathbf{b}$):</strong>
        $$3\\hat{x}_2 = ${qtb[1]} \\implies \\hat{x}_2 = ${x2}$$
        $$3\\hat{x}_1 + (${r12})(${x2}) = ${qtb[0]} \\implies 3\\hat{x}_1 = ${qtb[0] - r12 * x2} \\implies \\hat{x}_1 = ${x1}$$
      `
    };
  }

  function genInnerProductProblem() {
    if (Math.random() < 0.5) {
      const w1 = pick([1, 2, 3, 4]);
      const w2 = pick([1, 2, 3]);
      const a = pick([-2, -1, 1, 2]);
      const b = pick([-2, -1, 1, 2]);
      const c = pick([1, 2]);
      const d = pick([-2, -1, 1, 2]);
      const k = -(w1 * a * c + w2 * b * d);
      return {
        title: 'Weighted Inner Product Orthogonality',
        prompt: `Equip $\\mathbb{R}^3$ with the weighted inner product $\\langle \\mathbf{u}, \\mathbf{v} \\rangle = ${w1}u_1 v_1 + ${w2}u_2 v_2 + u_3 v_3$. Find the scalar $k$ so that $\\mathbf{u} = ${colVecTex([a, b, 1])}$ and $\\mathbf{v} = \\begin{bmatrix} ${c} \\\\ ${d} \\\\ k \\end{bmatrix}$ are orthogonal.`,
        fields: [{ id: 'k', label: 'k', ans: k }],
        solution: `
          Set $\\langle \\mathbf{u}, \\mathbf{v} \\rangle = 0$:
          $$${w1}(${a})(${c}) + ${w2}(${b})(${d}) + 1(1)(k) = 0$$
          $$${w1 * a * c} + (${w2 * b * d}) + k = 0 \\implies k = ${k}$$
        `
      };
    } else {
      const a = pick([-9, -6, -3, 3, 6, 9]);
      const k = -a / 3;
      return {
        title: 'Orthogonality in Function Space C[-1, 1]',
        prompt: `In $C[-1, 1]$ with inner product $\\langle f, g \\rangle = \\int_{-1}^1 f(t)g(t)\\,dt$, find the constant $k$ such that $p(t) = ${a}t^2 + k$ is orthogonal to $q(t) = 1$.`,
        fields: [{ id: 'k', label: 'k', ans: k }],
        solution: `
          Set $\\langle p, q \\rangle = \\int_{-1}^1 (${a}t^2 + k)\\cdot 1\\,dt = 0$:
          $$\\left[\\frac{${a}}{3}t^3 + kt\\right]_{-1}^{1} = \\frac{${2 * a}}{3} + 2k = ${2 * a / 3} + 2k = 0 \\implies k = ${k}$$
          <em>Connection:</em> Up to scaling, $p(t)$ is a multiple of the Legendre polynomial $P_2(t) = \\frac{1}{2}(3t^2 - 1)$!
        `
      };
    }
  }

  const GENERATORS = {
    proj: genProjectionProblem,
    gs: genGramSchmidtProblem,
    lsq: genLeastSquaresProblem,
    qr: genQRProblem,
    ips: genInnerProductProblem
  };

  // =========================================================================
  // 16-QUESTION MOCK EXAM BANK (5 SKILL AXES)
  // =========================================================================
  const AXES = [
    { key: 'norm', label: 'Norm & Orthog.', mod: 'Modules 1–2' },
    { key: 'proj', label: 'Proj. & Subspaces', mod: 'Modules 2–3' },
    { key: 'gsqr', label: 'Gram–Schmidt & QR', mod: 'Module 4' },
    { key: 'lsq',  label: 'Least Squares', mod: 'Modules 5–6' },
    { key: 'ips',  label: 'IPS & Fourier', mod: 'Module 7' }
  ];

  const EXAM_BANK = [
    {
      axis: 'norm',
      q: 'For two nonzero vectors $\\mathbf{u}, \\mathbf{v} \\in \\mathbb{R}^n$, when does equality $|\\mathbf{u}\\cdot\\mathbf{v}| = \\|\\mathbf{u}\\|\\|\\mathbf{v}\\|$ hold in the Cauchy–Schwarz inequality?',
      opts: [
        'If and only if $\\mathbf{u}$ and $\\mathbf{v}$ are orthogonal.',
        'If and only if $\\mathbf{u}$ and $\\mathbf{v}$ are linearly dependent (collinear).',
        'If and only if $\\|\\mathbf{u}\\| = \\|\\mathbf{v}\\| = 1$.',
        'Always, whenever $n \\ge 2$.'
      ],
      ans: 1,
      exp: 'Since $|\\mathbf{u}\\cdot\\mathbf{v}| = \\|\\mathbf{u}\\|\\|\\mathbf{v}\\|\\,|\\cos\\theta|$, equality holds iff $|\\cos\\theta| = 1$ ($\\theta = 0$ or $\\pi$), meaning one vector is a scalar multiple of the other.'
    },
    {
      axis: 'norm',
      q: 'Suppose $\\mathbf{u}, \\mathbf{v} \\in \\mathbb{R}^4$ satisfy $\\|\\mathbf{u}\\| = 3$, $\\|\\mathbf{v}\\| = 4$, and $\\|\\mathbf{u} + \\mathbf{v}\\| = 5$. What is $\\mathbf{u}\\cdot\\mathbf{v}$?',
      opts: ['$0$', '$5$', '$12$', '$-12$'],
      ans: 0,
      exp: '$\\|\\mathbf{u}+\\mathbf{v}\\|^2 = \\|\\mathbf{u}\\|^2 + 2(\\mathbf{u}\\cdot\\mathbf{v}) + \\|\\mathbf{v}\\|^2 \\implies 25 = 9 + 2(\\mathbf{u}\\cdot\\mathbf{v}) + 16 \\implies \\mathbf{u}\\cdot\\mathbf{v} = 0$ (Pythagorean Theorem!).'
    },
    {
      axis: 'norm',
      q: 'Let $U$ be an $m \\times n$ matrix ($m > n$) with <em>orthonormal columns</em>. Which statement is ALWAYS true?',
      opts: [
        '$U U^T = I_m$',
        '$U^T U = I_n$ and $\\|U\\mathbf{x}\\| = \\|\\mathbf{x}\\|$ for all $\\mathbf{x}\\in\\mathbb{R}^n$',
        '$U$ is invertible with $U^{-1} = U^T$',
        '$\\operatorname{Col} U = \\mathbb{R}^m$'
      ],
      ans: 1,
      exp: 'Because $m > n$, $U$ is rectangular (not invertible) and $U U^T$ is the rank-$n$ projection matrix onto $\\operatorname{Col} U \\ne I_m$. However, $U^T U = I_n$ ALWAYS holds and preserves lengths: $\\|U\\mathbf{x}\\|^2 = \\mathbf{x}^T U^T U \\mathbf{x} = \\|\\mathbf{x}\\|^2$.'
    },
    {
      axis: 'proj',
      q: 'If $A$ is a $5 \\times 8$ matrix with $\\operatorname{rank}(A) = 3$, what is the dimension of $(\\operatorname{Col} A)^\\perp$?',
      opts: ['$2$', '$3$', '$5$', '$8$'],
      ans: 0,
      exp: '$\\operatorname{Col} A$ is a subspace of $\\mathbb{R}^5$ with $\\dim(\\operatorname{Col} A) = \\operatorname{rank}(A) = 3$. Therefore $\\dim(\\operatorname{Col} A)^\\perp = \\dim(\\operatorname{Nul} A^T) = 5 - 3 = 2$.'
    },
    {
      axis: 'proj',
      q: 'Which fundamental subspace equality holds for EVERY $m \\times n$ matrix $A$?',
      opts: [
        '$(\\operatorname{Col} A)^\\perp = \\operatorname{Nul} A$',
        '$(\\operatorname{Row} A)^\\perp = \\operatorname{Nul} A^T$',
        '$(\\operatorname{Row} A)^\\perp = \\operatorname{Nul} A$',
        '$(\\operatorname{Nul} A)^\\perp = \\operatorname{Col} A$'
      ],
      ans: 2,
      exp: 'A vector $\\mathbf{x} \\in \\operatorname{Nul} A \\iff A\\mathbf{x} = \\mathbf{0} \\iff$ every row of $A$ is orthogonal to $\\mathbf{x} \\iff \\mathbf{x} \\in (\\operatorname{Row} A)^\\perp$.'
    },
    {
      axis: 'proj',
      q: 'Let $P$ be the orthogonal projection matrix onto a nontrivial subspace $W \\subset \\mathbb{R}^n$. Which pair of algebraic properties characterizes $P$?',
      opts: [
        '$P^2 = I$ and $P^T = -P$',
        '$P^2 = P$ and $P^T = P$',
        '$P^T P = I$ and $\\det P = 1$',
        '$P^2 = 0$ and $P^T = P$'
      ],
      ans: 1,
      exp: 'Projecting twice changes nothing ($P^2 = P$, idempotent), and orthogonal projections are always symmetric ($P^T = P$, e.g. $P = U U^T$).'
    },
    {
      axis: 'proj',
      q: 'If $\\hat{\\mathbf{y}} = \\operatorname{proj}_W(\\mathbf{y})$ and $\\mathbf{z} = \\mathbf{y} - \\hat{\\mathbf{y}}$, what is the minimum distance from $\\mathbf{y}$ to the subspace $W$?',
      opts: ['$\\|\\hat{\\mathbf{y}}\\|$', '$\\|\\mathbf{y}\\| - \\|\\hat{\\mathbf{y}}\\|$', '$\\|\\mathbf{z}\\|$', '$\\|\\mathbf{y}\\|^2 - \\|\\hat{\\mathbf{y}}\\|^2$'],
      ans: 2,
      exp: 'By the Best Approximation Theorem, $\\hat{\\mathbf{y}}$ is the unique closest point in $W$ to $\\mathbf{y}$, so the minimum distance is $\\|\\mathbf{y} - \\hat{\\mathbf{y}}\\| = \\|\\mathbf{z}\\|$.'
    },
    {
      axis: 'gsqr',
      q: 'When Gram–Schmidt converts a basis $\\{\\mathbf{x}_1, \\mathbf{x}_2, \\mathbf{x}_3\\}$ into an orthogonal basis $\\{\\mathbf{v}_1, \\mathbf{v}_2, \\mathbf{v}_3\\}$, which statement is ALWAYS true?',
      opts: [
        '$\\mathbf{v}_3$ is parallel to $\\mathbf{x}_3$.',
        '$\\operatorname{span}\\{\\mathbf{v}_1, \\mathbf{v}_2\\} = \\operatorname{span}\\{\\mathbf{x}_1, \\mathbf{x}_2\\}$.',
        'Every $\\mathbf{v}_i$ automatically has unit length $\\|\\mathbf{v}_i\\| = 1$.',
        'Reordering $\\{\\mathbf{x}_1, \\mathbf{x}_2, \\mathbf{x}_3\\}$ produces the exact same vectors $\\mathbf{v}_i$.'
      ],
      ans: 1,
      exp: 'At every stage $k$, Gram–Schmidt preserves the nested span: $\\operatorname{span}\\{\\mathbf{v}_1,\\dots,\\mathbf{v}_k\\} = \\operatorname{span}\\{\\mathbf{x}_1,\\dots,\\mathbf{x}_k\\}$. Normalization is an optional extra step, and reordering changes the resulting basis.'
    },
    {
      axis: 'gsqr',
      q: 'If $A$ is a $6 \\times 3$ matrix with linearly independent columns and $A = QR$ is its reduced $QR$ factorization, what are the dimensions and properties of $Q$ and $R$?',
      opts: [
        '$Q$ is $6 \\times 3$ (orthonormal columns); $R$ is $3 \\times 3$ upper triangular and invertible.',
        '$Q$ is $3 \\times 3$ orthogonal; $R$ is $6 \\times 3$ lower triangular.',
        '$Q$ is $6 \\times 6$ orthogonal; $R$ is $3 \\times 3$ diagonal.',
        '$Q$ is $6 \\times 3$ invertible; $R$ is $3 \\times 6$ upper triangular.'
      ],
      ans: 0,
      exp: 'In reduced $QR$, $Q$ has the same shape as $A$ ($6 \\times 3$) with orthonormal columns ($Q^T Q = I_3$), and $R = Q^T A$ is $3 \\times 3$ upper triangular with strictly positive diagonal entries (hence invertible).'
    },
    {
      axis: 'gsqr',
      q: 'In scientific computing, why is solving $R\\hat{\\mathbf{x}} = Q^T\\mathbf{b}$ preferred over directly solving the normal equations $A^T A\\hat{\\mathbf{x}} = A^T\\mathbf{b}$?',
      opts: [
        'Normal equations only work when $A$ is square.',
        'Forming $A^T A$ squares the condition number ($\\kappa(A^T A) = \\kappa(A)^2$), amplifying floating-point roundoff error.',
        '$R\\hat{\\mathbf{x}} = Q^T\\mathbf{b}$ gives a different vector $\\hat{\\mathbf{x}}$ with smaller residual.',
        '$A^T A$ is never symmetric.'
      ],
      ans: 1,
      exp: 'In exact arithmetic both yield the exact same $\\hat{\\mathbf{x}}$, but forming $A^T A$ squares the condition number, causing severe precision loss when columns of $A$ are nearly collinear.'
    },
    {
      axis: 'lsq',
      q: 'The least-squares problem $A\\mathbf{x} = \\mathbf{b}$ has a <strong>unique</strong> least-squares solution $\\hat{\\mathbf{x}}$ if and only if:',
      opts: [
        '$\\mathbf{b} \\in \\operatorname{Col} A$.',
        '$A$ is a square invertible matrix.',
        'The columns of $A$ are linearly independent (equivalently, $A^T A$ is invertible).',
        'The rows of $A$ are linearly independent.'
      ],
      ans: 2,
      exp: '$A^T A$ is invertible $\\iff \\operatorname{Nul}(A^T A) = \\operatorname{Nul}(A) = \\{\\mathbf{0}\\} \\iff$ the columns of $A$ are linearly independent.'
    },
    {
      axis: 'lsq',
      q: 'Geometrically, why do the normal equations $A^T A\\hat{\\mathbf{x}} = A^T\\mathbf{b}$ hold for a least-squares solution $\\hat{\\mathbf{x}}$?',
      opts: [
        'Because $A\\hat{\\mathbf{x}} = \\mathbf{b}$ is always consistent.',
        'Because $A^T(\\mathbf{b} - A\\hat{\\mathbf{x}}) = \\mathbf{0}$, meaning the residual $\\mathbf{b} - A\\hat{\\mathbf{x}}$ is orthogonal to $\\operatorname{Col} A$.',
        'Because $\\hat{\\mathbf{x}}$ must be orthogonal to $\\mathbf{b}$.',
        'Because $A^T A$ is always an orthogonal matrix.'
      ],
      ans: 1,
      exp: 'Minimizing $\\|\\mathbf{b} - A\\mathbf{x}\\|$ requires $A\\hat{\\mathbf{x}} = \\operatorname{proj}_{\\operatorname{Col} A}(\\mathbf{b})$, so $\\mathbf{b} - A\\hat{\\mathbf{x}} \\in (\\operatorname{Col} A)^\\perp = \\operatorname{Nul}(A^T)$, i.e. $A^T(\\mathbf{b} - A\\hat{\\mathbf{x}}) = \\mathbf{0}$.'
    },
    {
      axis: 'lsq',
      q: 'When fitting a quadratic model $y = \\beta_0 + \\beta_1 x + \\beta_2 x^2$ to $10$ data points $(x_i, y_i)$ by least squares, what is the size of the design matrix $X$?',
      opts: ['$3 \\times 10$', '$10 \\times 2$', '$10 \\times 3$', '$3 \\times 3$'],
      ans: 2,
      exp: 'Each of the $10$ observations contributes one row $[1,\\; x_i,\\; x_i^2]$, giving a $10 \\times 3$ design matrix $X$ and parameter vector $\\boldsymbol{\\beta} = (\\beta_0, \\beta_1, \\beta_2)^T$.'
    },
    {
      axis: 'ips',
      q: 'For which values of $w_1, w_2 \\in \\mathbb{R}$ does $\\langle \\mathbf{u}, \\mathbf{v} \\rangle = w_1 u_1 v_1 + w_2 u_2 v_2$ define a valid inner product on $\\mathbb{R}^2$?',
      opts: [
        'Any nonzero $w_1, w_2$',
        '$w_1 + w_2 > 0$',
        '$w_1 \\ge 0$ and $w_2 \\ge 0$',
        '$w_1 > 0$ and $w_2 > 0$'
      ],
      ans: 3,
      exp: 'Positivity requires $\\langle \\mathbf{u}, \\mathbf{u} \\rangle = w_1 u_1^2 + w_2 u_2^2 > 0$ for every nonzero $\\mathbf{u}$. Testing $\\mathbf{e}_1$ and $\\mathbf{e}_2$ shows both $w_1 > 0$ and $w_2 > 0$ are strictly required.'
    },
    {
      axis: 'ips',
      q: 'In $C[-1, 1]$ with $\\langle f, g \\rangle = \\int_{-1}^1 f(t)g(t)\\,dt$, which pair of functions is orthogonal?',
      opts: [
        '$f(t) = 1$ and $g(t) = t^2$',
        '$f(t) = t$ and $g(t) = t^2$',
        '$f(t) = 1$ and $g(t) = e^t$',
        '$f(t) = t$ and $g(t) = t^3$'
      ],
      ans: 1,
      exp: '$\\int_{-1}^1 t \\cdot t^2\\,dt = \\int_{-1}^1 t^3\\,dt = 0$ because $t^3$ is an odd function integrated over the symmetric interval $[-1, 1]$. All other pairs have strictly positive integrals.'
    },
    {
      axis: 'ips',
      q: 'In $C[-\\pi, \\pi]$, what geometric role does the $N$-th order Fourier partial sum $S_N(t)$ play?',
      opts: [
        'It matches the first $N$ derivatives of $f(t)$ at $t = 0$ (like a Taylor series).',
        'It interpolates $f(t)$ exactly at $N$ equally spaced points.',
        'It is the orthogonal projection of $f$ onto $\\operatorname{span}\\{1, \\cos t, \\sin t, \\dots, \\cos Nt, \\sin Nt\\}$, minimizing $\\int_{-\\pi}^\\pi |f(t) - g(t)|^2\\,dt$.',
        'It always has zero overshoot near jump discontinuities.'
      ],
      ans: 2,
      exp: 'Because the trigonometric basis is orthogonal under $\\langle f, g\\rangle = \\int_{-\\pi}^\\pi f(t)g(t)\\,dt$, $S_N(t)$ is literally the orthogonal projection onto the subspace $W_N$, minimizing the $L^2$ error.'
    }
  ];

  // =========================================================================
  // MODULE 8 REGISTRATION
  // =========================================================================
  AMA1751.registerModule(8, {
    render(meta) {
      const slideLabel = (meta && meta.slides) ? meta.slides : 'All Slides · Topic 6 Synthesis';
      return `
        <style>
          .m8-controls {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            margin-bottom: 12px;
            align-items: center;
          }
          .m8-chip {
            padding: 6px 12px;
            border-radius: 999px;
            border: 1px solid var(--border);
            background: #F8FAFC;
            font-size: 0.82rem;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.15s;
          }
          .m8-chip:hover, .m8-chip.active {
            background: var(--polyu-red-soft);
            border-color: #E8B4C4;
            color: var(--polyu-red-dark);
          }
          .m8-btn {
            padding: 8px 14px;
            border-radius: 8px;
            border: 1px solid #CBD5E1;
            background: #FFFFFF;
            font-weight: 700;
            font-size: 0.86rem;
            cursor: pointer;
            transition: all 0.15s;
          }
          .m8-btn:hover { background: #F1F5F9; }
          .m8-btn.primary {
            background: #0F766E;
            border-color: #0F766E;
            color: #FFFFFF;
          }
          .m8-btn.primary:hover { background: #115E59; }
          .m8-btn.accent {
            background: var(--polyu-red-dark);
            border-color: var(--polyu-red-dark);
            color: #FFFFFF;
          }
          .m8-atlas-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 12px;
          }
          .m8-atlas-card {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 12px 14px;
            font-size: 0.86rem;
            line-height: 1.55;
          }
          .m8-atlas-tag {
            display: inline-block;
            font-size: 0.7rem;
            font-weight: 800;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: #0F766E;
            background: #CCFBF1;
            padding: 2px 8px;
            border-radius: 999px;
            margin-bottom: 5px;
          }
          .m8-panel {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 14px 16px;
            font-size: 0.9rem;
            line-height: 1.6;
          }
          .m8-input-row {
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            align-items: center;
            margin: 12px 0;
          }
          .m8-field {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            background: #FFFFFF;
            border: 1px solid #CBD5E1;
            border-radius: 8px;
            padding: 4px 8px;
          }
          .m8-field label {
            font-weight: 700;
            font-size: 0.84rem;
            color: #334155;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          }
          .m8-field input {
            width: 66px;
            border: none;
            outline: none;
            font-size: 0.92rem;
            font-weight: 700;
            text-align: center;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          }
          .m8-badge {
            display: inline-block;
            padding: 5px 11px;
            border-radius: 8px;
            font-size: 0.84rem;
            font-weight: 700;
          }
          .m8-badge.ok { background: #DCFCE7; color: #166534; border: 1px solid #86EFAC; }
          .m8-badge.bad { background: #FEE2E2; color: #991B1B; border: 1px solid #FCA5A5; }
          .m8-badge.info { background: #E0F2FE; color: #075985; border: 1px solid #7DD3FC; }
          .m8-sol-box {
            margin-top: 12px;
            padding: 12px 14px;
            background: #FFFFFF;
            border-left: 4px solid #0F766E;
            border-radius: 0 8px 8px 0;
            border-top: 1px solid var(--border);
            border-right: 1px solid var(--border);
            border-bottom: 1px solid var(--border);
            font-size: 0.88rem;
          }
          .m8-exam-q {
            background: #F8FAFC;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 12px 14px;
            margin-bottom: 10px;
          }
          .m8-exam-opt {
            display: block;
            width: 100%;
            text-align: left;
            padding: 8px 12px;
            margin-top: 6px;
            border: 1px solid #CBD5E1;
            border-radius: 8px;
            background: #FFFFFF;
            cursor: pointer;
            font-size: 0.88rem;
            transition: all 0.14s;
          }
          .m8-exam-opt:hover { background: #F1F5F9; }
          .m8-exam-opt.selected { border-color: #1D4ED8; background: #EFF6FF; font-weight: 600; }
          .m8-exam-opt.correct { border-color: #16A34A; background: #DCFCE7; color: #166534; font-weight: 700; }
          .m8-exam-opt.wrong { border-color: #DC2626; background: #FEE2E2; color: #991B1B; }
          .m8-radar-wrap {
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 8px;
          }
          .m8-radar-canvas {
            display: block;
            width: 100%;
            height: 290px;
          }
          .m8-mat-grid {
            display: inline-grid;
            gap: 6px;
            padding: 8px 10px;
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 8px;
          }
          .m8-mat-grid input {
            width: 58px;
            padding: 5px;
            border: 1px solid #CBD5E1;
            border-radius: 6px;
            text-align: center;
            font-weight: 700;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          }
          .m8-stat-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(120px, 1fr));
            gap: 8px;
            margin: 10px 0;
          }
          .m8-stat-card {
            background: #FFFFFF;
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 8px 10px;
            text-align: center;
          }
          .m8-stat-lbl {
            font-size: 0.7rem;
            font-weight: 700;
            text-transform: uppercase;
            color: var(--text-muted);
          }
          .m8-stat-val {
            font-size: 0.96rem;
            font-weight: 800;
            color: #0F172A;
            margin-top: 2px;
            font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          }
        </style>

        <!-- HEADER -->
        <div class="module-header">
          <span class="slide-ref">Part 8 · ${slideLabel}</span>
          <h2>8. Topic 6 Mastery Arena, Procedural Drill &amp; Exam Simulator</h2>
          <p>Synthesize everything from <strong>Modules 1–7</strong>: review the unified Formula Atlas, train on infinite exact-integer exam problems, test yourself with a timed 10-question diagnostic mock exam, and verify your own homework via the universal $A = QR$ / Least-Squares Workbench.</p>
        </div>

        <!-- SECTION A: FILTERABLE TOPIC 6 FORMULA ATLAS -->
        <div class="card">
          <h3>Topic 6 Unified Formula Atlas</h3>
          <div class="m8-controls" id="m8AtlasFilters">
            <button class="m8-chip active" data-af="all">All Topics (6)</button>
            <button class="m8-chip" data-af="norm">Norm &amp; Orthogonality</button>
            <button class="m8-chip" data-af="proj">Projections &amp; Subspaces</button>
            <button class="m8-chip" data-af="gsqr">Gram–Schmidt &amp; QR</button>
            <button class="m8-chip" data-af="lsq">Least Squares</button>
            <button class="m8-chip" data-af="ips">Inner Product Spaces &amp; Fourier</button>
          </div>
          <div class="m8-atlas-grid" id="m8AtlasGrid"></div>
        </div>

        <!-- SECTION B: INTERACTIVE LAB 1 — INFINITE PROCEDURAL EXAM PROBLEM GENERATOR -->
        <div class="card">
          <h3>Interactive Lab 1: Infinite Exam Problem Drill (Exact-Integer Generator)</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            Every click generates a brand-new, algebraically clean problem with exact integer answers. Enter integers, decimals, or fractions (e.g. <code>-3/2</code>) and verify your steps.
          </p>

          <div class="m8-controls">
            <strong style="font-size:0.84rem;">Archetype:</strong>
            <button class="m8-chip active" data-gen="proj">Orthogonal Projection</button>
            <button class="m8-chip" data-gen="gs">Gram–Schmidt</button>
            <button class="m8-chip" data-gen="lsq">Normal Equations</button>
            <button class="m8-chip" data-gen="qr">QR Solver</button>
            <button class="m8-chip" data-gen="ips">Weighted / Integral IPS</button>
            <span class="m8-badge info" id="m8StreakBadge" style="margin-left:auto;">Streak: 0 🔥</span>
          </div>

          <div class="m8-panel">
            <div style="display:flex; justify-content:space-between; align-items:center; flex-wrap:wrap; gap:8px;">
              <strong id="m8ProbTitle" style="color:#0F172A; font-size:0.96rem;"></strong>
              <button class="m8-btn" id="m8NewProbBtn">↻ New Fresh Problem</button>
            </div>
            <div id="m8ProbPrompt" style="margin-top:8px;"></div>
            <div class="m8-input-row" id="m8ProbInputs"></div>
            <div style="display:flex; flex-wrap:wrap; gap:8px; align-items:center;">
              <button class="m8-btn primary" id="m8CheckProbBtn">Check Answer</button>
              <button class="m8-btn" id="m8ShowSolBtn">Show Step-by-Step Solution</button>
              <span id="m8ProbFeedback" style="font-weight:700; font-size:0.88rem;"></span>
            </div>
            <div class="m8-sol-box" id="m8ProbSol" style="display:none;"></div>
          </div>
        </div>

        <!-- SECTION C: INTERACTIVE LAB 2 — ADAPTIVE MOCK EXAM SIMULATOR + RADAR CHART -->
        <div class="card">
          <h3>Interactive Lab 2: Topic 6 Diagnostic Mock Exam (10 Questions + Mastery Radar)</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            Draws a balanced 10-question exam (2 questions per competency axis) from the 16-question pitfall bank. Switch between untimed <strong>Practice Mode</strong> and <strong>15-Minute Timed Exam</strong>.
          </p>

          <div class="m8-controls">
            <button class="m8-btn primary" id="m8StartPracticeBtn">Start Practice Exam (Untimed)</button>
            <button class="m8-btn accent" id="m8StartTimedBtn">⏱ Start Timed Exam (15:00)</button>
            <span class="m8-badge info" id="m8TimerBadge" style="display:none;">Time Left: 15:00</span>
            <span class="m8-badge ok" id="m8ScoreBadge" style="display:none;"></span>
          </div>

          <div class="grid-2" style="align-items:start;">
            <div>
              <div id="m8ExamContainer"></div>
              <div id="m8ExamActions" style="margin-top:10px; display:none;">
                <button class="m8-btn primary" id="m8SubmitExamBtn">Grade Exam &amp; Update Radar</button>
              </div>
            </div>

            <div class="m8-panel">
              <div style="font-weight:700; margin-bottom:6px;">5-Axis Topic 6 Mastery Radar:</div>
              <div class="m8-radar-wrap">
                <canvas id="m8RadarCanvas" class="m8-radar-canvas"></canvas>
              </div>
              <div id="m8RadarSummary" style="margin-top:10px; font-size:0.86rem; color:#334155;">
                Complete and grade an exam to populate your personalized mastery profile and study prescription.
              </div>
            </div>
          </div>
        </div>

        <!-- SECTION D: INTERACTIVE LAB 3 — UNIVERSAL TOPIC 6 HOMEWORK VERIFIER WORKBENCH -->
        <div class="card">
          <h3>Interactive Lab 3: Universal Topic 6 Pipeline Workbench ($A \\in \\mathbb{R}^{3\\times 2},\\; \\mathbf{b} \\in \\mathbb{R}^3$)</h3>
          <p style="margin-top:0; color:var(--text-muted);">
            Verify your own tutorial or homework problems! Edit any entry of $A$ and $\\mathbf{b}$ to compute <strong>Gram–Schmidt</strong>, <strong>Reduced $QR$</strong>, <strong>Normal Equations</strong>, and the <strong>Orthogonal Decomposition $\\mathbf{b} = \\hat{\\mathbf{b}} + \\mathbf{e}$</strong> simultaneously.
          </p>

          <div class="m8-controls">
            <strong style="font-size:0.84rem;">Presets:</strong>
            <button class="m8-chip active" data-wb="classic">Classic Full-Rank</button>
            <button class="m8-chip" data-wb="linefit">Line Fit (3 Data Points)</button>
            <button class="m8-chip" data-wb="orth">Already Orthogonal Columns</button>
            <button class="m8-chip" data-wb="singular">Rank-Deficient Trap (Singular $A^T A$!)</button>
          </div>

          <div style="display:flex; flex-wrap:wrap; gap:18px; align-items:center; margin-bottom:14px;">
            <div>
              <div style="font-size:0.82rem; font-weight:700; margin-bottom:4px;">Matrix $A$ ($3\\times 2$):</div>
              <div class="m8-mat-grid" style="grid-template-columns: repeat(2, 1fr);" id="m8WbA"></div>
            </div>
            <div>
              <div style="font-size:0.82rem; font-weight:700; margin-bottom:4px;">Vector $\\mathbf{b}$ ($3\\times 1$):</div>
              <div class="m8-mat-grid" style="grid-template-columns: 1fr;" id="m8WbB"></div>
            </div>
          </div>

          <div class="m8-panel" id="m8WbOutput"></div>

          <details class="code-details">
            <summary>Python (NumPy) Companion Code — One-Shot Topic 6 Verification Pipeline</summary>
            <pre><code>import numpy as np

A = np.array([[1., 0.],
              [1., 1.],
              [0., 1.]])
b = np.array([2., 1., 2.])

# 1) Reduced QR Factorization (A = Q @ R)
Q, R = np.linalg.qr(A, mode='reduced')
print("Q^T Q =\n", np.round(Q.T @ Q, 6))
print("R =\n", np.round(R, 4))

# 2) Least-Squares Solution via QR vs Normal Equations
x_qr = np.linalg.solve(R, Q.T @ b)
x_ne = np.linalg.solve(A.T @ A, A.T @ b)
print("x_hat (QR) =", np.round(x_qr, 4), "| (Normal Eq) =", np.round(x_ne, 4))

# 3) Orthogonal Decomposition b = b_hat + e & Pythagorean Check
b_hat = A @ x_qr
e = b - b_hat
print("b_hat . e =", np.round(np.dot(b_hat, e), 10))
print("||b_hat||^2 + ||e||^2 =", np.dot(b_hat, b_hat) + np.dot(e, e), "| ||b||^2 =", np.dot(b, b))</code></pre>
          </details>
        </div>
      `;
    },

    mount(container, { renderMath }) {
      // =====================================================================
      // SECTION A: FILTERABLE FORMULA ATLAS
      // =====================================================================
      const atlasGrid = container.querySelector('#m8AtlasGrid');
      function renderAtlas(filter = 'all') {
        const list = filter === 'all' ? FORMULAS : FORMULAS.filter(f => f.cat === filter);
        atlasGrid.innerHTML = list.map(f => `
          <div class="m8-atlas-card">
            <div class="m8-atlas-tag">${f.tag}</div>
            <div style="font-weight:800; color:#0F172A; margin-bottom:4px;">${f.title}</div>
            <div>${f.body}</div>
          </div>
        `).join('');
        renderMath(atlasGrid);
      }
      container.querySelectorAll('[data-af]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-af]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          renderAtlas(btn.dataset.af);
        });
      });
      renderAtlas('all');

      // =====================================================================
      // SECTION B: LAB 1 — PROCEDURAL EXAM PROBLEM ENGINE
      // =====================================================================
      let activeGen = 'proj';
      let currentProb = null;
      let streak = 0;
      let probSolved = false;

      function loadNewProblem() {
        currentProb = GENERATORS[activeGen]();
        probSolved = false;
        container.querySelector('#m8ProbTitle').textContent = currentProb.title;
        const promptEl = container.querySelector('#m8ProbPrompt');
        promptEl.innerHTML = currentProb.prompt;
        const inpEl = container.querySelector('#m8ProbInputs');
        inpEl.innerHTML = currentProb.fields.map(f => `
          <span class="m8-field">
            <label for="m8Inp_${f.id}">${f.label} =</label>
            <input type="text" id="m8Inp_${f.id}" placeholder="?" autocomplete="off" />
          </span>
        `).join('');
        const fb = container.querySelector('#m8ProbFeedback');
        fb.textContent = '';
        const sol = container.querySelector('#m8ProbSol');
        sol.style.display = 'none';
        sol.innerHTML = currentProb.solution;
        renderMath(promptEl);
      }

      container.querySelectorAll('[data-gen]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-gen]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          activeGen = btn.dataset.gen;
          loadNewProblem();
        });
      });

      container.querySelector('#m8NewProbBtn').addEventListener('click', loadNewProblem);

      container.querySelector('#m8CheckProbBtn').addEventListener('click', () => {
        const fb = container.querySelector('#m8ProbFeedback');
        let allOk = true;
        let anyBlank = false;
        currentProb.fields.forEach(f => {
          const el = container.querySelector(`#m8Inp_${f.id}`);
          const val = parseNum(el.value);
          if (!Number.isFinite(val)) {
            anyBlank = true;
            allOk = false;
            el.style.color = '#DC2626';
          } else if (Math.abs(val - f.ans) < 1e-4) {
            el.style.color = '#166534';
          } else {
            allOk = false;
            el.style.color = '#DC2626';
          }
        });
        if (anyBlank) {
          fb.style.color = '#B45309';
          fb.textContent = 'Please enter valid numbers (or fractions like -3/2) in all boxes.';
          return;
        }
        if (allOk) {
          if (!probSolved) {
            streak += 1;
            probSolved = true;
            container.querySelector('#m8StreakBadge').textContent = `Streak: ${streak} 🔥`;
          }
          fb.style.color = '#166534';
          fb.textContent = '✓ Exact match! Click "New Fresh Problem" to keep your streak going.';
          const sol = container.querySelector('#m8ProbSol');
          sol.style.display = 'block';
          renderMath(sol);
        } else {
          streak = 0;
          container.querySelector('#m8StreakBadge').textContent = 'Streak: 0 🔥';
          fb.style.color = '#991B1B';
          fb.textContent = '✗ Not quite — red entries differ from the exact answer. Try again or inspect the solution!';
        }
      });

      container.querySelector('#m8ShowSolBtn').addEventListener('click', () => {
        const sol = container.querySelector('#m8ProbSol');
        sol.style.display = sol.style.display === 'none' ? 'block' : 'none';
        if (sol.style.display === 'block') renderMath(sol);
      });

      loadNewProblem();

      // =====================================================================
      // SECTION C: LAB 2 — ADAPTIVE MOCK EXAM SIMULATOR + RADAR CHART
      // =====================================================================
      let examQs = [];
      let examAnswers = {};
      let examGraded = false;
      let timerId = null;
      let secondsLeft = 900;
      let radarScores = [0, 0, 0, 0, 0];

      const radarCanvas = container.querySelector('#m8RadarCanvas');
      const radarCtx = radarCanvas.getContext('2d');

      function drawRadar(scores, hasData = false) {
        if (!radarCanvas.isConnected) return;
        const dpr = window.devicePixelRatio || 1;
        const rect = radarCanvas.getBoundingClientRect();
        const w = rect.width || 360, h = rect.height || 290;
        radarCanvas.width = Math.round(w * dpr);
        radarCanvas.height = Math.round(h * dpr);
        radarCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
        radarCtx.clearRect(0, 0, w, h);

        const cx = w / 2, cy = h / 2 + 8;
        const R = Math.min(w, h) * 0.34;
        const N = AXES.length;
        const angleAt = i => -Math.PI / 2 + (2 * Math.PI * i) / N;

        // Concentric pentagons (25%, 50%, 75%, 100%)
        [0.25, 0.5, 0.75, 1.0].forEach(frac => {
          radarCtx.beginPath();
          for (let i = 0; i < N; i++) {
            const a = angleAt(i);
            const x = cx + R * frac * Math.cos(a);
            const y = cy + R * frac * Math.sin(a);
            if (i === 0) radarCtx.moveTo(x, y); else radarCtx.lineTo(x, y);
          }
          radarCtx.closePath();
          radarCtx.strokeStyle = frac === 1 ? '#CBD5E1' : '#E2E8F0';
          radarCtx.lineWidth = frac === 1 ? 1.5 : 1;
          radarCtx.stroke();
        });

        // Radial spokes & labels
        radarCtx.font = 'bold 11px sans-serif';
        for (let i = 0; i < N; i++) {
          const a = angleAt(i);
          const ex = cx + R * Math.cos(a), ey = cy + R * Math.sin(a);
          radarCtx.beginPath();
          radarCtx.moveTo(cx, cy); radarCtx.lineTo(ex, ey);
          radarCtx.strokeStyle = '#E2E8F0'; radarCtx.lineWidth = 1; radarCtx.stroke();

          const lx = cx + (R + 24) * Math.cos(a);
          const ly = cy + (R + 18) * Math.sin(a);
          radarCtx.fillStyle = '#334155';
          radarCtx.textAlign = Math.cos(a) > 0.25 ? 'left' : (Math.cos(a) < -0.25 ? 'right' : 'center');
          radarCtx.textBaseline = 'middle';
          const pctStr = hasData ? ` (${Math.round(scores[i] * 100)}%)` : '';
          radarCtx.fillText(`${AXES[i].label}${pctStr}`, lx, ly);
        }

        // Data polygon
        if (hasData) {
          radarCtx.beginPath();
          for (let i = 0; i < N; i++) {
            const a = angleAt(i);
            const rVal = Math.max(0.08, scores[i]) * R;
            const x = cx + rVal * Math.cos(a), y = cy + rVal * Math.sin(a);
            if (i === 0) radarCtx.moveTo(x, y); else radarCtx.lineTo(x, y);
          }
          radarCtx.closePath();
          radarCtx.fillStyle = 'rgba(15, 118, 110, 0.24)';
          radarCtx.fill();
          radarCtx.strokeStyle = '#0F766E';
          radarCtx.lineWidth = 2.6;
          radarCtx.stroke();

          for (let i = 0; i < N; i++) {
            const a = angleAt(i);
            const rVal = Math.max(0.08, scores[i]) * R;
            radarCtx.beginPath();
            radarCtx.arc(cx + rVal * Math.cos(a), cy + rVal * Math.sin(a), 4.5, 0, Math.PI * 2);
            radarCtx.fillStyle = scores[i] < 0.5 ? '#DC2626' : '#0F766E';
            radarCtx.fill();
          }
        }
      }

      function stopTimer() {
        if (timerId) { clearInterval(timerId); timerId = null; }
      }

      function formatTime(sec) {
        const m = Math.floor(sec / 60);
        const s = sec % 60;
        return `${String(m).padStart(2, '0')}:${String(s).padStart(2, '0')}`;
      }

      function buildBalancedExam() {
        // Pick 2 questions from each of the 5 axes (10 total)
        const picked = [];
        AXES.forEach(ax => {
          const pool = shuffle(EXAM_BANK.filter(item => item.axis === ax.key));
          picked.push(...pool.slice(0, 2));
        });
        return shuffle(picked);
      }

      function renderExam() {
        const box = container.querySelector('#m8ExamContainer');
        box.innerHTML = examQs.map((item, idx) => {
          const chosen = examAnswers[idx];
          const optsHtml = item.opts.map((opt, oi) => {
            let cls = 'm8-exam-opt';
            if (!examGraded && chosen === oi) cls += ' selected';
            if (examGraded) {
              if (oi === item.ans) cls += ' correct';
              else if (chosen === oi && oi !== item.ans) cls += ' wrong';
            }
            const letter = String.fromCharCode(65 + oi);
            return `<button class="${cls}" data-qidx="${idx}" data-oidx="${oi}" ${examGraded ? 'disabled' : ''}><strong>${letter}.</strong> ${opt}</button>`;
          }).join('');
          const expHtml = examGraded ? `
            <div style="margin-top:8px; font-size:0.83rem; color:${chosen === item.ans ? '#166534' : '#991B1B'};">
              <strong>${chosen === item.ans ? '✓ Correct!' : '✗ Review:'}</strong> ${item.exp}
            </div>
          ` : '';
          return `
            <div class="m8-exam-q">
              <div style="font-size:0.89rem;"><strong>Q${idx + 1}.</strong> ${item.q}</div>
              ${optsHtml}
              ${expHtml}
            </div>
          `;
        }).join('');

        if (!examGraded) {
          box.querySelectorAll('.m8-exam-opt').forEach(btn => {
            btn.addEventListener('click', () => {
              const qi = parseInt(btn.dataset.qidx, 10);
              const oi = parseInt(btn.dataset.oidx, 10);
              examAnswers[qi] = oi;
              // Update button visual state cleanly without full re-render flicker
              const parent = btn.closest('.m8-exam-q');
              parent.querySelectorAll('.m8-exam-opt').forEach(b => b.classList.remove('selected'));
              btn.classList.add('selected');
            });
          });
        }
        renderMath(box);
      }

      function gradeExam() {
        if (examGraded || examQs.length === 0) return;
        examGraded = true;
        stopTimer();

        let totalCorrect = 0;
        const axisHits = {}, axisTotal = {};
        AXES.forEach(a => { axisHits[a.key] = 0; axisTotal[a.key] = 0; });

        examQs.forEach((q, idx) => {
          axisTotal[q.axis] += 1;
          if (examAnswers[idx] === q.ans) {
            totalCorrect += 1;
            axisHits[q.axis] += 1;
          }
        });

        radarScores = AXES.map(a => axisTotal[a.key] ? axisHits[a.key] / axisTotal[a.key] : 0);
        drawRadar(radarScores, true);
        renderExam();

        const scoreBadge = container.querySelector('#m8ScoreBadge');
        scoreBadge.style.display = 'inline-block';
        scoreBadge.className = totalCorrect >= 7 ? 'm8-badge ok' : 'm8-badge bad';
        scoreBadge.textContent = `Score: ${totalCorrect} / ${examQs.length} (${Math.round((totalCorrect / examQs.length) * 100)}%)`;

        container.querySelector('#m8ExamActions').style.display = 'none';

        // Diagnostic prescription
        const minScore = Math.min(...radarScores);
        const weakest = AXES.filter((_, i) => radarScores[i] === minScore);
        const sumEl = container.querySelector('#m8RadarSummary');
        if (totalCorrect === examQs.length) {
          sumEl.innerHTML = `<div class="m8-badge ok">🏆 Topic 6 Mastered (100%)!</div><p style="margin:8px 0 0;">Flawless across all 5 axes. Challenge yourself with Timed Mode or fresh procedural drills above.</p>`;
        } else {
          const weakList = weakest.map(w => `<strong>${w.label}</strong> (${w.mod})`).join(', ');
          sumEl.innerHTML = `
            <div><strong>Diagnostic Prescription:</strong></div>
            <p style="margin:4px 0 0;">Focus your revision on ${weakList}. Review the worked explanations highlighted in red on the left, then run 3 fresh problems in Lab 1!</p>
          `;
        }
      }

      function startExam(timed) {
        stopTimer();
        examQs = buildBalancedExam();
        examAnswers = {};
        examGraded = false;
        container.querySelector('#m8ScoreBadge').style.display = 'none';
        container.querySelector('#m8ExamActions').style.display = 'block';
        const tBadge = container.querySelector('#m8TimerBadge');

        if (timed) {
          secondsLeft = 900;
          tBadge.style.display = 'inline-block';
          tBadge.textContent = `Time Left: ${formatTime(secondsLeft)}`;
          timerId = setInterval(() => {
            if (!radarCanvas.isConnected) { stopTimer(); return; }
            secondsLeft -= 1;
            tBadge.textContent = `Time Left: ${formatTime(Math.max(0, secondsLeft))}`;
            if (secondsLeft <= 0) gradeExam();
          }, 1000);
        } else {
          tBadge.style.display = 'none';
        }
        renderExam();
      }

      container.querySelector('#m8StartPracticeBtn').addEventListener('click', () => startExam(false));
      container.querySelector('#m8StartTimedBtn').addEventListener('click', () => startExam(true));
      container.querySelector('#m8SubmitExamBtn').addEventListener('click', gradeExam);

      // Start with an initial ready-to-take Practice Exam & empty radar
      drawRadar(radarScores, false);
      startExam(false);

      // =====================================================================
      // SECTION D: LAB 3 — UNIVERSAL TOPIC 6 PIPELINE WORKBENCH
      // =====================================================================
      const wbPresets = {
        classic:  { A: [[1, 0], [1, 1], [0, 1]],   b: [2, 1, 2] },
        linefit:  { A: [[1, 0], [1, 1], [1, 2]],   b: [1, 2, 4] },
        orth:     { A: [[1, 1], [1, -1], [0, 0]],  b: [3, 1, 2] },
        singular: { A: [[1, 2], [2, 4], [-1, -2]], b: [1, 0, 3] }
      };

      let wbA = wbPresets.classic.A.map(r => r.slice());
      let wbB = wbPresets.classic.b.slice();

      function buildWbInputs() {
        const aWrap = container.querySelector('#m8WbA');
        const bWrap = container.querySelector('#m8WbB');
        aWrap.innerHTML = '';
        bWrap.innerHTML = '';
        for (let r = 0; r < 3; r++) {
          for (let c = 0; c < 2; c++) {
            const inp = document.createElement('input');
            inp.type = 'number'; inp.step = 'any'; inp.value = wbA[r][c];
            inp.addEventListener('input', () => {
              wbA[r][c] = parseFloat(inp.value) || 0;
              computeWb();
            });
            aWrap.appendChild(inp);
          }
          const binp = document.createElement('input');
          binp.type = 'number'; binp.step = 'any'; binp.value = wbB[r];
          binp.addEventListener('input', () => {
            wbB[r] = parseFloat(binp.value) || 0;
            computeWb();
          });
          bWrap.appendChild(binp);
        }
      }

      function dot3(u, v) { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; }

      function computeWb() {
        const out = container.querySelector('#m8WbOutput');
        const a1 = [wbA[0][0], wbA[1][0], wbA[2][0]];
        const a2 = [wbA[0][1], wbA[1][1], wbA[2][1]];
        const b = wbB.slice();

        const ata00 = dot3(a1, a1), ata01 = dot3(a1, a2), ata11 = dot3(a2, a2);
        const detATA = ata00 * ata11 - ata01 * ata01;
        const atb = [dot3(a1, b), dot3(a2, b)];

        if (ata00 < 1e-10 || detATA < 1e-9) {
          out.innerHTML = `
            <div class="m8-badge bad">✗ Columns of $A$ are Linearly Dependent ($\\det(A^T A) = ${clean(detATA, 4)} \\approx 0$)</div>
            <p style="margin:8px 0 0; font-size:0.86rem; color:#7F1D1D;">
              <strong>Exam Trap Detected!</strong> Because $\\mathbf{a}_1$ and $\\mathbf{a}_2$ are collinear (or zero), Gram–Schmidt produces $\\mathbf{v}_2 = \\mathbf{0}$, reduced $QR$ does not have an invertible $R$, and the normal equations $A^T A \\hat{\\mathbf{x}} = A^T \\mathbf{b}$ have <strong>infinitely many</strong> least-squares solutions rather than a unique $\\hat{\\mathbf{x}}$.
            </p>
          `;
          renderMath(out);
          return;
        }

        // 1) Gram-Schmidt & QR
        const v1 = a1.slice();
        const n1 = Math.sqrt(dot3(v1, v1));
        const q1 = v1.map(x => x / n1);
        const r11 = n1;
        const r12 = dot3(q1, a2);
        const v2 = a2.map((x, i) => x - r12 * q1[i]);
        const n2 = Math.sqrt(dot3(v2, v2));
        const q2 = v2.map(x => x / n2);
        const r22 = n2;

        // 2) Solve via QR: R x_hat = Q^T b
        const qtb1 = dot3(q1, b), qtb2 = dot3(q2, b);
        const x2 = qtb2 / r22;
        const x1 = (qtb1 - r12 * x2) / r11;

        // 3) Orthogonal decomposition b = b_hat + e
        const bHat = [0, 1, 2].map(i => x1 * a1[i] + x2 * a2[i]);
        const err = [0, 1, 2].map(i => b[i] - bHat[i]);
        const nB2 = dot3(b, b), nBH2 = dot3(bHat, bHat), nE2 = dot3(err, err);
        const orthoCheck = dot3(bHat, err);

        const Qmat = [[q1[0], q2[0]], [q1[1], q2[1]], [q1[2], q2[2]]];
        const Rmat = [[r11, r12], [0, r22]];
        const ATAmat = [[ata00, ata01], [ata01, ata11]];

        out.innerHTML = `
          <div class="m8-badge ok">✓ Full Column Rank ($\\det(A^T A) = ${clean(detATA, 3)} > 0$) — Unique Least-Squares Solution Exists</div>
          <div class="grid-2" style="margin-top:10px; gap:12px;">
            <div style="background:#FFFFFF; border:1px solid var(--border); border-radius:8px; padding:10px 12px;">
              <strong>1. Gram–Schmidt &amp; Reduced $QR$:</strong>
              <div style="overflow-x:auto; margin-top:4px;">
                $$\\mathbf{v}_1 = ${colVecTex(v1)}, \\quad \\mathbf{v}_2 = ${colVecTex(v2)}$$
                $$Q = ${matTex(Qmat)}, \\quad R = ${matTex(Rmat)}$$
              </div>
            </div>
            <div style="background:#FFFFFF; border:1px solid var(--border); border-radius:8px; padding:10px 12px;">
              <strong>2. Normal Equations &amp; Solution:</strong>
              <div style="overflow-x:auto; margin-top:4px;">
                $$A^T A = ${matTex(ATAmat)}, \\quad A^T \\mathbf{b} = ${colVecTex(atb)}$$
                $$\\hat{\\mathbf{x}} = ${colVecTex([x1, x2])}$$
              </div>
            </div>
          </div>

          <div style="background:#FFFFFF; border:1px solid var(--border); border-radius:8px; padding:10px 12px; margin-top:10px;">
            <strong>3. Orthogonal Decomposition $\\mathbf{b} = \\hat{\\mathbf{b}} + \\mathbf{e}$ ($\\hat{\\mathbf{b}} \\in \\operatorname{Col} A,\\; \\mathbf{e} \\in \\operatorname{Nul} A^T$):</strong>
            <div style="overflow-x:auto; margin-top:4px;">
              $$${colVecTex(b)} = \\underbrace{${colVecTex(bHat)}}_{\\hat{\\mathbf{b}} = A\\hat{\\mathbf{x}}} + \\underbrace{${colVecTex(err)}}_{\\mathbf{e} = \\mathbf{b} - A\\hat{\\mathbf{x}}}$$
            </div>
            <div class="m8-stat-grid">
              <div class="m8-stat-card"><div class="m8-stat-lbl">$\\hat{\\mathbf{b}} \\cdot \\mathbf{e}$</div><div class="m8-stat-val" style="color:#166534;">${clean(orthoCheck, 6)}</div></div>
              <div class="m8-stat-card"><div class="m8-stat-lbl">Error $\\|\\mathbf{e}\\|$</div><div class="m8-stat-val" style="color:#D97706;">${clean(Math.sqrt(nE2), 4)}</div></div>
              <div class="m8-stat-card"><div class="m8-stat-lbl">$\\|\\hat{\\mathbf{b}}\\|^2 + \\|\\mathbf{e}\\|^2$</div><div class="m8-stat-val">${clean(nBH2 + nE2, 3)}</div></div>
              <div class="m8-stat-card"><div class="m8-stat-lbl">$\\|\\mathbf{b}\\|^2$</div><div class="m8-stat-val">${clean(nB2, 3)}</div></div>
            </div>
          </div>
        `;
        renderMath(out);
      }

      container.querySelectorAll('[data-wb]').forEach(btn => {
        btn.addEventListener('click', () => {
          container.querySelectorAll('[data-wb]').forEach(b => b.classList.remove('active'));
          btn.classList.add('active');
          const p = wbPresets[btn.dataset.wb];
          wbA = p.A.map(r => r.slice());
          wbB = p.b.slice();
          buildWbInputs();
          computeWb();
        });
      });

      buildWbInputs();
      computeWb();

      window.addEventListener('resize', () => {
        drawRadar(radarScores, examGraded);
      });
    }
  });
})();
