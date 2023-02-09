"use strict";
(() => {
  var __defProp = Object.defineProperty;
  var __defProps = Object.defineProperties;
  var __getOwnPropDescs = Object.getOwnPropertyDescriptors;
  var __getOwnPropSymbols = Object.getOwnPropertySymbols;
  var __hasOwnProp = Object.prototype.hasOwnProperty;
  var __propIsEnum = Object.prototype.propertyIsEnumerable;
  var __defNormalProp = (obj, key, value) => key in obj ? __defProp(obj, key, { enumerable: true, configurable: true, writable: true, value }) : obj[key] = value;
  var __spreadValues = (a, b) => {
    for (var prop in b ||= {})
      if (__hasOwnProp.call(b, prop))
        __defNormalProp(a, prop, b[prop]);
    if (__getOwnPropSymbols)
      for (var prop of __getOwnPropSymbols(b)) {
        if (__propIsEnum.call(b, prop))
          __defNormalProp(a, prop, b[prop]);
      }
    return a;
  };
  var __spreadProps = (a, b) => __defProps(a, __getOwnPropDescs(b));

  // src/CanvasPen.ts
  var CanvasPen = class {
    constructor(ctx) {
      this.ctx = ctx;
    }
    moveTo(x, y) {
      this.ctx.moveTo(x, y);
      this.lastx = x;
      this.lasty = y;
    }
    arcTo(x, y, turn) {
      if (this.lastx == void 0 || this.lasty == void 0) {
        throw new Error("CanvasPen - arc/line with no current point");
      }
      const absTurn = Math.abs(turn);
      if (absTurn < 1e-5) {
        this.ctx.lineTo(x, y);
      } else {
        const sinHalf = Math.sin(absTurn * 0.5);
        const tanHalf = sinHalf / Math.cos(absTurn * 0.5);
        const yfac = turn <= 0 ? 1 : -1;
        const dx = x - this.lastx;
        const dy = y - this.lasty;
        const px = this.lastx + (dx - dy * tanHalf * yfac) * 0.5;
        const py = this.lasty + (dy + dx * tanHalf * yfac) * 0.5;
        const r = Math.sqrt(dx * dx + dy * dy) * 0.5 / sinHalf;
        this.ctx.arcTo(px, py, x, y, r);
      }
      this.lastx = x;
      this.lasty = y;
    }
  };

  // src/biarc.ts
  function biArcSplit(pt0, pt1) {
    const [x0, y0, tx0, ty0] = pt0;
    const [x1, y1, tx1, ty1] = pt1;
    const dx = x1 - x0;
    const dy = y1 - y0;
    let dmag = Math.sqrt(dx * dx + dy * dy);
    const txm = dx / dmag;
    const tym = dy / dmag;
    const f0 = ty1 * dx - tx1 * dy;
    const f1 = tx0 * dy - ty0 * dx;
    const k = dmag / ((tx0 * txm + ty0 * tym + 1) * f0 + (tx1 * txm + ty1 * tym + 1) * f1);
    const a0 = k * f0;
    const xm = x0 + (tx0 + txm) * a0;
    const ym = y0 + (ty0 + tym) * a0;
    return [xm, ym, txm, tym];
  }
  function biArcApproximate(samples, tolerance) {
    const best = [[0, 0, -1]];
    let nextScanStart = 0;
    for (let pos = 1; pos < samples.length; ++pos) {
      const thisScanStart = nextScanStart;
      nextScanStart = pos - 1;
      let bestCount = best[pos - 1][0] + 1;
      let bestMaxError = best[pos - 1][1];
      let bestPred = pos - 1;
      for (let testPos = thisScanStart; testPos < pos - 1; ++testPos) {
        const error = getApproxError(samples, testPos, pos);
        if (error > tolerance) {
          continue;
        }
        if (testPos < nextScanStart) {
          nextScanStart = testPos;
        }
        const count = best[testPos][0] + 1;
        const maxError = Math.max(best[testPos][1], error);
        if (count > bestCount) {
          break;
        }
        if (count < bestCount || maxError < bestMaxError) {
          bestCount = count;
          bestMaxError = maxError;
          bestPred = testPos;
        }
        if (error < maxError) {
          break;
        }
      }
      best.push([bestCount, bestMaxError, bestPred]);
    }
    const ret = [];
    for (let i = best.length - 1; i >= 0; i = best[i][2]) {
      ret.push(samples[i]);
    }
    return ret.reverse();
  }
  function getApproxError(samples, fromPos, toPos) {
    const pt0 = samples[fromPos];
    const pt1 = samples[toPos];
    const ptm = biArcSplit(pt0, pt1);
    let maxError = 0;
    for (let i = fromPos + 1; i < toPos; ++i) {
      maxError = Math.max(maxError, getBiArcPointError(pt0, ptm, pt1, samples[i]));
    }
    return maxError;
  }
  function getBiArcPointError(pt0, ptm, pt1, sample) {
    const [x0, y0] = pt0;
    const [x1, y1] = pt1;
    const [xm, ym] = ptm;
    const [xs, ys] = sample;
    const dx = x1 - x0;
    const dy = y1 - y0;
    if ((xs - x0) * dx + (ys - y0) * dy > (xm - x0) * dx + (ym - y0) * dy) {
      return getArcPointError(ptm, pt1, sample);
    } else {
      return getArcPointError(pt0, ptm, sample);
    }
  }
  function getArcPointError(pt0, pt1, sample) {
    const [x0, y0, tx0, ty0] = pt0;
    const [x1, y1, tx1, ty1] = pt1;
    const [xs, ys] = sample;
    const dxa = x1 - x0;
    const dya = y1 - y0;
    const den = dxa * (ty1 - ty0) + dya * (tx0 - tx1);
    const da2 = dxa * dxa + dya * dya;
    if (da2 > Math.abs(den) * 1e8) {
      const xh = (x0 + x1) * 0.5;
      const yh = (y0 + y1) * 0.5;
      return Math.abs((xs - xh) * (ty0 + ty1) - (ys - yh) * (tx0 + tx1)) * 0.5;
    }
    const r = da2 / den;
    const xc = (x0 + x1 - r * (ty0 + ty1)) * 0.5;
    const yc = (y0 + y1 + r * (tx0 + tx1)) * 0.5;
    const xcs = xs - xc;
    const ycs = ys - yc;
    const rs = Math.sqrt(xcs * xcs + ycs * ycs);
    return Math.abs(rs - Math.abs(r));
  }
  function drawBiArc(pen, pt0, pt1) {
    const [, , tx0, ty0] = pt0;
    const [x1, y1, tx1, ty1] = pt1;
    const [xm, ym, txm, tym] = biArcSplit(pt0, pt1);
    const rot0 = Math.asin(tx0 * tym - ty0 * txm);
    const rot1 = Math.asin(txm * ty1 - tym * tx1);
    pen.arcTo(xm, ym, rot0);
    pen.arcTo(x1, y1, rot1);
  }

  // src/floatBinarySearch.ts
  function searchForFloat(lo, hi, isLowerBound) {
    if (!(lo < hi)) {
      return [lo, hi];
    }
    [lo, hi] = getLinearRange(lo, hi, isLowerBound);
    for (; ; ) {
      const testVal = lo + (hi - lo) * 0.5;
      if (testVal <= lo || testVal >= hi) {
        break;
      }
      if (isLowerBound(testVal)) {
        lo = testVal;
      } else {
        hi = testVal;
      }
    }
    return [lo, hi];
  }
  function getLinearRange(low, high, isLowerBound) {
    let negRange;
    if (low < 0) {
      if (high > 0) {
        if (isLowerBound(0)) {
          return scaleRange(0, high, 0.25, isLowerBound);
        } else {
          const isNegLowerBound = (n) => !isLowerBound(-n);
          negRange = scaleRange(0, -low, 0.25, isNegLowerBound);
        }
      } else {
        const isNegLowerBound = (n) => !isLowerBound(-n);
        negRange = scaleRange(-high, -low, 0.25, isNegLowerBound);
      }
    } else {
      return scaleRange(low, high, 0.25, isLowerBound);
    }
    low = -negRange[1];
    negRange[1] = -negRange[0];
    negRange[0] = low;
    return negRange;
  }
  function scaleRange(low, high, minScale, isLowerBound) {
    if (!(minScale > 0 && low < high * minScale)) {
      return [low, high];
    }
    const range = scaleRange(low, high, minScale * minScale, isLowerBound);
    [low, high] = range;
    const test = high * minScale;
    if (test > low && test < high) {
      if (isLowerBound(test)) {
        range[0] = test;
      } else {
        range[1] = test;
      }
    }
    return range;
  }

  // src/CircleCut.ts
  var CircleCut = class {
    constructor(a0, x0, y0, a1, x1, y1, arcTol) {
      this.arcTolerance = arcTol;
      if (a1 < a0) {
        [a0, a1] = [a1, a0];
        [x0, x1] = [x1, x0];
        [y0, y1] = [y1, y0];
      }
      this.x0 = x0;
      this.x1 = x1;
      this.y0 = y0;
      this.y1 = y1;
      this.a0 = a0;
      this.a1 = a1;
      const dx = x1 - x0;
      const dy = y1 - y0;
      const da = a1 - a0;
      this.DNUM = x0 * dy - y0 * dx;
      this.DA = dx * dx + dy * dy;
      this.DB = 2 * (x0 * dx + y0 * dy);
      this.mag02 = x0 * x0 + y0 * y0;
      const DC = this.mag02 - this.DNUM / da;
      const disc = this.DB * this.DB - 4 * this.DA * DC;
      this.tmidpoint = -this.DB / (2 * this.DA);
      this.thetaMidpoint = this.getThetaForT(this.tmidpoint);
      if (disc > 0) {
        const root = Math.sqrt(disc);
        const tstart = Math.max((-this.DB - root) / (2 * this.DA), 0);
        const tend = Math.min((-this.DB + root) / (2 * this.DA), 1);
        if (tstart < tend) {
          this.reversal = {
            tstart,
            tend,
            thetaStart: this.getThetaForT(tstart),
            thetaEnd: this.getThetaForT(tend)
          };
        }
      }
    }
    getDiscontinuityThetas(minTheta, maxTheta) {
      const ret = [];
      if (!this.reversal) {
        return ret;
      }
      for (const theta of [this.reversal.thetaStart, this.reversal.thetaEnd]) {
        if (theta > minTheta && theta < maxTheta) {
          ret.push(theta);
        }
      }
      return ret;
    }
    drawSegment(pen, thetaFrom, thetaTo, doInitialMove) {
      let tstart;
      let tend;
      if (this.reversal) {
        if (thetaTo > this.reversal.thetaEnd) {
          if (thetaFrom < this.reversal.thetaEnd) {
            this.drawSegment(pen, thetaFrom, this.reversal.thetaEnd, doInitialMove);
            thetaFrom = this.reversal.thetaEnd;
            doInitialMove = false;
          }
          tstart = this.getTForThetaInBase(thetaFrom);
          tend = this.getTForThetaInBase(thetaTo);
        } else if (thetaTo > this.reversal.thetaStart) {
          if (thetaFrom < this.reversal.thetaStart) {
            this.drawSegment(pen, thetaFrom, this.reversal.thetaStart, doInitialMove);
            thetaFrom = this.reversal.thetaStart;
            doInitialMove = false;
          }
          tstart = this.getTForThetaInReversal(thetaFrom);
          tend = this.getTForThetaInReversal(thetaTo);
        } else {
          tstart = this.getTForThetaInBase(thetaFrom);
          tend = this.getTForThetaInBase(thetaTo);
        }
      } else if (thetaTo > this.thetaMidpoint) {
        if (thetaFrom < this.thetaMidpoint) {
          this.drawSegment(pen, thetaFrom, this.thetaMidpoint, doInitialMove);
          thetaFrom = this.thetaMidpoint;
          doInitialMove = false;
        }
        tstart = this.getTForThetaInBase(thetaFrom);
        tend = this.getTForThetaInBase(thetaTo);
      } else {
        tstart = this.getTForThetaInBase(thetaFrom);
        tend = this.getTForThetaInBase(thetaTo);
      }
      const tangentSign = Math.sign(tend - tstart);
      const pt0 = this.getPointAndTangentForT(tstart, tangentSign);
      const pt1 = this.getPointAndTangentForT(tend, tangentSign);
      if (doInitialMove) {
        pen.moveTo(pt0[0], pt0[1]);
      } else {
        pen.arcTo(pt0[0], pt0[1], 0);
      }
      let samples = [pt0];
      this.getApproxSamples(samples, tstart, pt0, tend, pt1);
      samples = biArcApproximate(samples, this.arcTolerance);
      for (let i = 1; i < samples.length; ++i) {
        drawBiArc(pen, samples[i - 1], samples[i]);
      }
    }
    getApproxSamples(out, tstart, ptstart, tend, ptend) {
      if (tend == tstart) {
        return;
      }
      const tmid = tstart + (tend - tstart) * 0.5;
      if (tmid == tstart || tmid == tend) {
        out.push(ptend);
        return;
      }
      const tangentSign = Math.sign(tend - tstart);
      const ptmid = this.getPointAndTangentForT(tmid, tangentSign);
      const [x0, y0] = ptstart;
      const [x1, y1] = ptend;
      const [xm, ym] = ptmid;
      const dx = x1 - x0;
      const dy = y1 - y0;
      let dev = (xm - x0) * dy - (ym - y0) * dx;
      dev = dev * dev / (dx * dx + dy * dy);
      if (dev > this.arcTolerance * this.arcTolerance * 0.25) {
        this.getApproxSamples(out, tstart, ptstart, tmid, ptmid);
        this.getApproxSamples(out, tmid, ptmid, tend, ptend);
      } else {
        out.push(ptend);
      }
    }
    getR(theta) {
      return this.getRForT(this.getTForTheta(theta));
    }
    getThetaForT(t) {
      const invt = 1 - t;
      const y = this.y1 + (this.y0 - this.y1) * invt;
      const x = this.x1 + (this.x0 - this.x1) * invt;
      const a = this.a1 + (this.a0 - this.a1) * invt;
      return Math.atan2(y, x) - a;
    }
    getRForT(t) {
      const invt = 1 - t;
      const y = this.y1 + (this.y0 - this.y1) * invt;
      const x = this.x1 + (this.x0 - this.x1) * invt;
      return Math.sqrt(x * x + y * y);
    }
    getPointAndTangentForT(t, tangentSign) {
      const invt = 1 - t;
      const dpxdt = this.x1 - this.x0;
      const dpydt = this.y1 - this.y0;
      const dadt = this.a1 - this.a0;
      const py = this.y1 - dpydt * invt;
      const px = this.x1 - dpxdt * invt;
      const rotc = Math.cos(dadt * invt - this.a1);
      const rots = Math.sin(dadt * invt - this.a1);
      const x = px * rotc - py * rots;
      const y = py * rotc + px * rots;
      let dxdt = dpxdt * rotc - dpydt * rots + y * dadt;
      let dydt = dpydt * rotc + dpxdt * rots - x * dadt;
      let mag2 = dxdt * dxdt + dydt * dydt;
      if (mag2 < 1e-16) {
        dxdt = x;
        dydt = y;
        mag2 = dxdt * dxdt + dydt * dydt;
        if (t < this.tmidpoint) {
          tangentSign = -tangentSign;
        }
      }
      const normfac = tangentSign / Math.sqrt(mag2);
      return [x, y, dxdt * normfac, dydt * normfac];
    }
    getTForTheta(theta) {
      if (this.reversal && theta >= this.reversal.thetaStart && theta <= this.reversal.thetaEnd) {
        return this.getTForThetaInReversal(theta);
      }
      return this.getTForThetaInBase(theta);
    }
    getTForThetaInBase(theta) {
      let tmin = 0;
      let tmax = 1;
      if (this.reversal) {
        if (theta < this.thetaMidpoint) {
          tmin = this.reversal.tend;
        } else {
          tmax = this.reversal.tstart;
        }
      }
      [tmin, tmax] = searchForFloat(tmin, tmax, (t) => this.getThetaForT(t) >= theta);
      const thetaMin = this.getThetaForT(tmin);
      const thetaMax = this.getThetaForT(tmax);
      let interp = 0;
      if (thetaMax != thetaMin) {
        interp = Math.max((theta - thetaMin) / (thetaMax - thetaMin), 0);
        interp = Math.min(interp, 1);
      }
      return tmin + (tmax - tmin) * interp;
    }
    getTForThetaInReversal(theta) {
      if (!this.reversal) {
        return this.tmidpoint;
      }
      let tmin = this.reversal.tstart;
      let tmax = this.reversal.tend;
      [tmin, tmax] = searchForFloat(tmin, tmax, (t) => this.getThetaForT(t) <= theta);
      const thetaMin = this.getThetaForT(tmin);
      const thetaMax = this.getThetaForT(tmax);
      let interp = 0.5;
      if (thetaMax != thetaMin) {
        interp = Math.max((theta - thetaMin) / (thetaMax - thetaMin), 0);
        interp = Math.min(interp, 1);
      }
      return tmin + (tmax - tmin) * interp;
    }
  };

  // src/ConstantRadiusCut.ts
  var ConstantRadiusCut = class {
    constructor(r) {
      this.r = r;
    }
    getDiscontinuityThetas(_minTheta, _maxTheta) {
      return [];
    }
    drawSegment(pen, thetaFrom, thetaTo, doInitialMove) {
      if (thetaFrom - thetaTo > Math.PI * 0.6) {
        const mid = thetaFrom + (thetaTo - thetaFrom) * 0.5;
        this.drawSegment(pen, thetaFrom, mid, doInitialMove);
        this.drawSegment(pen, mid, thetaTo, false);
        return;
      }
      const sx = Math.cos(thetaFrom) * this.r;
      const sy = Math.sin(thetaFrom) * this.r;
      const ex = Math.cos(thetaTo) * this.r;
      const ey = Math.sin(thetaTo) * this.r;
      if (doInitialMove) {
        pen.moveTo(sx, sy);
      } else {
        pen.arcTo(sx, sy, 0);
      }
      pen.arcTo(ex, ey, thetaTo - thetaFrom);
    }
    getR(theta) {
      return this.r;
    }
  };

  // src/pq.ts
  function pqpush(q, val) {
    let pos = q.length;
    q.push(val);
    while (pos > 0) {
      let parpos = (pos - 1) / 2 | 0;
      if (q[parpos][0] <= q[pos][0]) {
        break;
      }
      let t = q[parpos];
      q[parpos] = q[pos];
      q[pos] = t;
      pos = parpos;
    }
  }
  function pqpop(q) {
    if (!q.length) {
      return void 0;
    }
    if (q.length < 2) {
      return q.shift();
    }
    let ret = q[0];
    q[0] = q.pop();
    let pos = 0;
    let leftpos = 1;
    while (leftpos < q.length) {
      let minpos = pos;
      if (q[minpos][0] > q[leftpos][0]) {
        minpos = leftpos;
      }
      let rightpos = leftpos + 1;
      if (rightpos < q.length && q[minpos][0] > q[rightpos][0]) {
        minpos = rightpos;
      }
      if (minpos == pos) {
        break;
      }
      let t = q[pos];
      q[pos] = q[minpos];
      q[minpos] = t;
      pos = minpos;
      leftpos = minpos * 2 + 1;
    }
    return ret;
  }

  // src/pathSampling.ts
  var BOTTOM_TOLERANCE = 1e-5;
  function normalizePolarCutPath(segments, dadt) {
    const pq = [];
    for (let cut of segments) {
      let [sa, ea] = cut;
      if (ea == sa) {
        continue;
      }
      if (ea < sa) {
        sa = cut[1];
        ea = cut[0];
      }
      const curve = cut[2];
      let rot = cut[3];
      if (sa < -0.5 || sa >= 0.5) {
        const flr = Math.floor(sa + 0.5);
        sa -= flr;
        ea -= flr;
        rot -= flr;
      }
      while (ea > 0.5) {
        pqpush(pq, [sa, 0.5, curve, rot]);
        sa = -0.5;
        ea -= 1;
        rot -= 1;
      }
      pqpush(pq, [sa, ea, curve, rot]);
    }
    if (!pq.length) {
      return [];
    }
    const sampleAs = [];
    {
      const eventPointsSet = /* @__PURE__ */ new Set();
      for (const [sa, ea, curve, rot] of pq) {
        eventPointsSet.add(sa);
        eventPointsSet.add(ea);
        for (const a of curve.getDiscontinuityThetas((sa - rot) * dadt, (ea - rot) * dadt)) {
          eventPointsSet.add(a / dadt + rot);
        }
      }
      const eventAs = [...eventPointsSet].sort((a, b) => a - b);
      const avoidance = 1e-6;
      let rangeStart = -0.5;
      for (const a of eventAs) {
        if (a - avoidance >= rangeStart) {
          const rangeEnd = a - avoidance;
          if (rangeEnd - rangeStart < avoidance) {
            sampleAs.push(rangeStart + (rangeEnd - rangeStart) * 0.5);
          }
          const n = Math.floor((rangeEnd - rangeStart) / 1e-3) + 1;
          for (let i = 0; i <= n; ++i) {
            sampleAs.push(rangeStart + (rangeEnd - rangeStart) * i / n);
          }
        }
        rangeStart = a + avoidance;
      }
    }
    let prevCandidates = [];
    let prevCandidateStartSample = -0.5;
    let bottomCuts = [];
    let bottomCutStartSamples = [];
    let bottomCutEndSamples = [];
    let preva = -0.5;
    let activeCuts = [];
    for (const cura of sampleAs) {
      while (pq.length && pq[0][0] <= cura) {
        activeCuts.push(pqpop(pq));
      }
      let d = 0;
      for (let s = 0; s < activeCuts.length; ++s) {
        if (activeCuts[s][1] > cura) {
          activeCuts[d++] = activeCuts[s];
        }
      }
      activeCuts.length = d;
      const curMin = getMinR(activeCuts, dadt, cura);
      if (prevCandidates.length) {
        const survivors = curMin == null ? [] : getBottomSegs(prevCandidates, dadt, cura, curMin);
        if (!survivors.length) {
          bottomCuts.push(prevCandidates[0]);
          bottomCutStartSamples.push(prevCandidateStartSample);
          bottomCutEndSamples.push(preva);
        }
        prevCandidates = survivors;
        preva = cura;
      }
      if (!prevCandidates.length && curMin != null) {
        prevCandidates = getBottomSegs(activeCuts, dadt, cura, curMin);
        prevCandidateStartSample = cura;
      }
    }
    if (prevCandidates.length) {
      bottomCuts.push(prevCandidates[0]);
      bottomCutStartSamples.push(prevCandidateStartSample);
      bottomCutEndSamples.push(preva);
    }
    for (let i = 1; i < bottomCuts.length; ++i) {
      const locut = bottomCuts[i - 1];
      const hicut = bottomCuts[i];
      if (locut[1] > hicut[0]) {
        let loa = Math.max(bottomCutStartSamples[i - 1], hicut[0]);
        let hia = Math.min(bottomCutStartSamples[i], locut[1]);
        [loa, hia] = searchForFloat(loa, hia, (testa) => {
          const lor2 = locut[2].getR((testa - locut[3]) * dadt);
          const hir2 = hicut[2].getR((testa - hicut[3]) * dadt);
          return lor2 < hir2;
        });
        const lor = locut[2].getR((loa - locut[3]) * dadt);
        const hir = hicut[2].getR((hia - hicut[3]) * dadt);
        if (lor == hir) {
          hia = loa;
        }
        locut[1] = loa;
        hicut[0] = hia;
      }
    }
    return bottomCuts;
  }
  function getBottomSegs(cuts, dadt, a, bottom) {
    let ret = [];
    for (const cut of cuts) {
      if (cut[1] >= a) {
        const r = cut[2].getR((a - cut[3]) * dadt);
        if (r <= bottom + BOTTOM_TOLERANCE) {
          ret.push(cut);
        }
      }
    }
    return ret;
  }
  function getMinR(cuts, dadt, a) {
    let ret = null;
    for (const cut of cuts) {
      const r = cut[2].getR((a - cut[3]) * dadt);
      if (ret == null || r <= ret) {
        ret = r;
      }
    }
    return ret;
  }

  // src/XFormPen.ts
  var XFormPen = class {
    constructor(delegate) {
      this.delegate = delegate;
      this.rotDegrees = 0;
      this.flipYinFac = 1;
      this.scaleFactor = 1;
      this.xx = 1;
      this.xy = 0;
      this.tx = 0;
      this.ty = 0;
    }
    rotate(degrees) {
      this.rotDegrees += degrees * this.flipYinFac;
      this.rotDegrees -= Math.floor(this.rotDegrees / 360) * 360;
      this.xx = Math.cos(this.rotDegrees * Math.PI / 180);
      this.xy = Math.sin(this.rotDegrees * Math.PI / 180);
      const quarters = this.rotDegrees / 90;
      if (Math.floor(quarters) === quarters) {
        if ((quarters & 1) == 0) {
          this.xx = Math.sign(this.xx);
          this.xy = 0;
        } else {
          this.xx = 0;
          this.xy = Math.sign(this.xy);
        }
      }
      this.xx *= this.scaleFactor;
      this.xy *= this.scaleFactor;
      return this;
    }
    transformPoint(x, y) {
      return [
        this.tx + x * this.xx - y * this.flipYinFac * this.xy,
        this.ty + x * this.xy + y * this.flipYinFac * this.xx
      ];
    }
    translate(x, y) {
      [this.tx, this.ty] = this.transformPoint(x, y);
      return this;
    }
    scale(fac, flipY) {
      if (fac < 0) {
        this.rotate(180);
        this.scaleFactor *= -fac;
      } else {
        this.scaleFactor *= fac;
      }
      if (flipY) {
        this.flipYinFac = -this.flipYinFac;
      }
      this.xx = Math.cos(this.rotDegrees * Math.PI / 180) * this.scaleFactor;
      this.xy = Math.sin(this.rotDegrees * Math.PI / 180) * this.scaleFactor;
      return this;
    }
    copy() {
      return new XFormPen(this.delegate).translate(this.tx, this.ty).rotate(this.rotDegrees).scale(this.scaleFactor, this.flipYinFac < 0);
    }
    moveTo(x, y) {
      const [newX, newY] = this.transformPoint(x, y);
      this.delegate.moveTo(newX, newY);
    }
    arcTo(x, y, leftTurnRadians) {
      const [newX, newY] = this.transformPoint(x, y);
      this.delegate.arcTo(newX, newY, leftTurnRadians * this.flipYinFac);
    }
  };

  // src/GearCutter.ts
  var GearCutter = class {
    constructor(nTeeth, pitchRadius, faceTol, filletTol) {
      this.lastPointIsCut = false;
      this.pointCurves = /* @__PURE__ */ new Map();
      this.flatCurves = /* @__PURE__ */ new Map();
      this.path = [];
      this.nTeeth = nTeeth;
      this.pitchRadius = pitchRadius;
      this.dadTooth = Math.PI * 2 / nTeeth;
      this.dydTooth = this.dadTooth * pitchRadius;
      this.faceTol = faceTol;
      this.filletTol = filletTol;
    }
    drawToothPath(pen, doInitialMove) {
      const segments = normalizePolarCutPath(this.path, this.dadTooth);
      for (const seg of segments) {
        const [sa, ea, c, rot] = seg;
        const xpen = new XFormPen(pen).rotate(rot * 360 / this.nTeeth);
        c.drawSegment(xpen, (sa - rot) * this.dadTooth, (ea - rot) * this.dadTooth, doInitialMove);
        doInitialMove = false;
      }
    }
    moveTo(x, y) {
      if (x <= 0) {
        throw new Error("x <= 0 is not supported in GearCutter");
      }
      this.lastX = x;
      this.lastY = y;
      this.lastPointIsCut = false;
    }
    arcTo(x, y, turn) {
      if (x <= 0) {
        throw new Error("x <= 0 is not supported in GearCutter");
      }
      if (Math.abs(turn) > 1e-3) {
        throw new Error("Curved cutter paths are not supported in GearCutter");
      }
      if (this.lastX == void 0 || this.lastY == void 0) {
        throw new Error("Curve without current point sent to GearCutter");
      }
      let x0 = this.lastX;
      let y0 = this.lastY;
      this.lastX = x;
      this.lastY = y;
      if (!this.lastPointIsCut) {
        this.cutPoint(x0, y0);
      }
      this.lastPointIsCut = true;
      if (x0 == x && y0 == y) {
        return;
      }
      this.cutPoint(x, y);
      if (x0 == x) {
        this.cutFlat(x, y0, y);
        return;
      }
      const xp = this.pitchRadius;
      const y0p = (y - y0) * (xp - x0) / (x - x0) + y0;
      const tp = -y0p / this.dydTooth;
      let dirx = y0 - y;
      let diry = x - x0;
      let k = this.dydTooth * diry / (dirx * dirx + diry * diry);
      let dxdt = k * dirx;
      let dydt = k * diry;
      const t0 = (x0 - xp) / dxdt;
      const t1 = (x - xp) / dxdt;
      const curve = new CircleCut(
        (t0 + tp) * this.dadTooth,
        x0,
        t0 * dydt,
        (t1 + tp) * this.dadTooth,
        x,
        t1 * dydt,
        this.faceTol
      );
      const thetaA = curve.getThetaForT(0) / this.dadTooth;
      const thetaB = curve.getThetaForT(1) / this.dadTooth;
      this.path.push([Math.min(thetaA, thetaB), Math.max(thetaA, thetaB), curve, 0]);
    }
    cutPoint(x, y) {
      let curve = this.pointCurves.get(x);
      if (!curve) {
        curve = new CircleCut(
          -Math.PI,
          x,
          -Math.PI * this.pitchRadius,
          Math.PI,
          x,
          Math.PI * this.pitchRadius,
          this.filletTol
        );
        this.pointCurves.set(x, curve);
      }
      const rot = y / this.dydTooth;
      this.path.push([rot - Math.PI / this.dadTooth, rot + Math.PI / this.dadTooth, curve, rot]);
    }
    cutFlat(x, y0, y1) {
      let curve = this.flatCurves.get(x);
      if (!curve) {
        curve = new ConstantRadiusCut(x);
        this.flatCurves.set(x, curve);
      }
      if (y0 < y1) {
        this.path.push([y0 / this.dydTooth, y1 / this.dydTooth, curve, 0]);
      } else {
        this.path.push([y1 / this.dydTooth, y0 / this.dydTooth, curve, 0]);
      }
    }
  };

  // src/rack.ts
  function makeRack(props) {
    const {
      contactRatio: contactRatio2,
      pressureAngle: pressureAngle2,
      profileShift: profileShift2,
      balancePercent: balancePercent2,
      topClrPercent: topReliefPercent,
      botClrPercent: botReliefPercent,
      balanceAbsPercent
    } = props;
    const sinPA = Math.sin(pressureAngle2 * Math.PI / 180);
    const cosPA = Math.cos(pressureAngle2 * Math.PI / 180);
    const tanPA = sinPA / cosPA;
    let ah = contactRatio2 * sinPA * cosPA;
    let cy = profileShift2 / (100 * Math.PI);
    let miny = cy - ah / 2;
    let maxy = cy + ah / 2;
    let bkw = balanceAbsPercent / (200 * Math.PI);
    let freew = 0.5 - ah * tanPA;
    const cx = -0.25 - freew * (balancePercent2 - 50) / 100;
    maxy += topReliefPercent / (100 * Math.PI);
    miny -= botReliefPercent / (100 * Math.PI);
    const topx = (maxy - cy) * tanPA + cx;
    const botx = (miny - cy) * tanPA + cx;
    return (pen, doMove) => {
      if (doMove) {
        pen.moveTo(-1 - botx + bkw, miny);
      }
      pen.arcTo(botx - bkw, miny, 0);
      pen.arcTo(topx - bkw, maxy, 0);
      pen.arcTo(-topx + bkw, maxy, 0);
      pen.arcTo(-botx + bkw, miny, 0);
    };
  }

  // src/RecordingPen.ts
  var RecordingPen = class {
    constructor() {
      this.xs = [];
      this.ys = [];
      this.turns = [];
      this.path = (pen, doMove) => {
        const { xs, ys, turns } = this;
        let i = 0;
        if (!doMove) {
          while (i < turns.length && turns[i] == null) {
            ++i;
          }
        }
        for (; i < turns.length; ++i) {
          const x = xs[i];
          const y = ys[i];
          const a = turns[i];
          if (a == null) {
            pen.moveTo(x, y);
          } else {
            pen.arcTo(x, y, a);
          }
        }
      };
    }
    moveTo(x, y) {
      if (this.turns.length && this.turns[this.turns.length - 1] == null) {
        this.turns.pop();
        this.xs.pop();
        this.ys.pop();
      }
      this.xs.push(x);
      this.ys.push(y);
      this.turns.push(null);
    }
    arcTo(x, y, turn) {
      if (!this.turns.length) {
        throw new Error("arc without preceding move in RecordingPen");
      } else {
        const i = this.turns.length - 1;
        const dx = x - this.xs[i];
        const dy = y - this.ys[i];
        const mag2 = dx * dx + dy * dy;
        if (mag2 < 1e-14) {
          return;
        }
        if (mag2 < 1e-8) {
          turn = 0;
        }
      }
      this.xs.push(x);
      this.ys.push(y);
      this.turns.push(turn);
    }
    countSegments() {
      let ret = 0;
      for (const t of this.turns) {
        if (t != null) {
          ++ret;
        }
      }
      return ret;
    }
  };

  // src/svg.ts
  var SvgRecorder = class {
    constructor(props) {
      this.buf = [
        '<?xml version="1.0" encoding="UTF-8"?>\n',
        "<svg>"
      ];
      this.bounds = null;
      this.scale = props.scale || -96 / 25.4;
    }
    getSvgText() {
      const oldlen = this.buf.length;
      let [minx, miny, width, height] = this.bounds || [0, 0, 0, 0];
      minx = Math.floor(minx - 0.1);
      miny = Math.floor(miny - 0.1);
      width = Math.ceil(width + 0.1) - minx;
      height = Math.ceil(height + 0.1) - miny;
      this.buf[1] = `<svg viewBox="${minx} ${miny} ${width} ${height}" width="${width}" height="${height}" version="1.1" xmlns="http://www.w3.org/2000/svg">
`;
      this.buf.push("</svg>\n");
      const svg = this.buf.join("");
      this.buf.length = oldlen;
      return svg;
    }
    drawCircle(drawProps, cx, cy, radius) {
      const fill = drawProps.fill || "none";
      const stroke = drawProps.stroke || "none";
      const strokeWidth = (drawProps.strokeWidth || 1) * Math.abs(this.scale);
      cx *= Math.abs(this.scale);
      cy *= this.scale;
      radius *= Math.abs(this.scale);
      this.buf.push(`<g fill="${fill}" stroke-width="${strokeWidth}" stroke="${stroke}">
`);
      this.buf.push(`    <circle cx="${cx}" cy="${cy}" r="${radius}"/>
`);
      this.buf.push("</g>");
      this.mergeBounds([cx - radius, cy - radius, cx + radius, cy + radius]);
    }
    draw(drawProps, path) {
      var _a;
      const pen = new SvgPen(
        __spreadProps(__spreadValues({}, drawProps), {
          strokeWidth: ((_a = drawProps.strokeWidth) != null ? _a : 1) * Math.abs(this.scale)
        }),
        this.buf
      );
      const xfpen = new XFormPen(pen);
      xfpen.scale(Math.abs(this.scale), this.scale < 0);
      path(xfpen, true);
      this.mergeBounds(pen.finish());
    }
    mergeBounds(bounds) {
      if (!this.bounds) {
        this.bounds = bounds;
      } else if (bounds) {
        for (let i = 0; i < 2; ++i) {
          this.bounds[i] = Math.min(this.bounds[i], bounds[i]);
          this.bounds[i + 2] = Math.max(this.bounds[i + 2], bounds[i + 2]);
        }
      }
    }
  };
  var SvgPen = class {
    constructor(drawProps, buf) {
      this.havePoint = false;
      this.havePath = false;
      this.lastx = 0;
      this.lasty = 0;
      this.minx = 0;
      this.miny = 0;
      this.maxx = 0;
      this.maxy = 0;
      this.buf = buf;
      this.pendingProps = drawProps;
      this.closePaths = !!drawProps.fill;
    }
    finish() {
      if (this.pendingProps) {
        return null;
      }
      if (this.havePath) {
        this.buf.push(this.closePaths ? ' Z"/>\n' : '"/>\n');
        this.havePath = false;
      }
      this.buf.push("</g>");
      return [this.minx, this.miny, this.maxx, this.maxy];
    }
    moveTo(x, y) {
      if (this.havePath) {
        this.buf.push(this.closePaths ? ' Z"/>\n' : '"/>\n');
        this.havePath = false;
      }
      this.lastx = x;
      this.lasty = y;
      this.havePoint = true;
    }
    arcTo(x, y, turn) {
      var _a;
      if (!this.havePoint) {
        throw new Error("arcTo before moveTo in SVG with");
      }
      if (this.pendingProps) {
        const fill = this.pendingProps.fill || "none";
        const stroke = this.pendingProps.stroke || "none";
        const strokeWidth = (_a = this.pendingProps.strokeWidth) != null ? _a : 1;
        this.pendingProps = null;
        this.buf.push(`<g fill="${fill}" stroke-width="${strokeWidth}" stroke="${stroke}">
`);
        this.minx = Math.min(this.lastx, x);
        this.miny = Math.min(this.lasty, y);
        this.maxx = Math.max(this.lastx, x);
        this.maxy = Math.max(this.lasty, y);
      } else {
        this.minx = Math.min(this.minx, Math.min(this.lastx, x));
        this.maxx = Math.max(this.maxx, Math.max(this.lastx, x));
        this.miny = Math.min(this.miny, Math.min(this.lasty, y));
        this.maxy = Math.max(this.maxy, Math.max(this.lasty, y));
      }
      if (!this.havePath) {
        this.buf.push(`    <path d=" M ${this.lastx} ${this.lasty}`);
        this.havePath = true;
      }
      if (Math.abs(turn) < 1e-5) {
        this.buf.push(` L ${x} ${y}`);
      } else {
        const side = turn > 0 ? 1 : 0;
        const dx = x - this.lastx;
        const dy = y - this.lasty;
        const dist = Math.sqrt(dx * dx + dy * dy);
        const a = Math.abs(turn * 0.5);
        const sina = Math.sin(a);
        const r = dist * 0.5 / sina;
        const tanfac = (turn >= 0 ? 0.5 : -0.5) / Math.tan(a);
        const cx = x - dx * 0.5 - dy * tanfac;
        const cy = y - dy * 0.5 + dx * tanfac;
        if (cx > Math.min(this.lastx, x) && cx < Math.max(this.lastx, x)) {
          if (cy > (y + this.lasty) * 0.5) {
            this.miny = Math.min(this.miny, cy - r);
          } else {
            this.maxy = Math.max(this.maxy, cy + r);
          }
        }
        if (cy > Math.min(this.lasty, y) && cy < Math.max(this.lasty, y)) {
          if (cx > (x + this.lastx) * 0.5) {
            this.minx = Math.min(this.minx, cx - r);
          } else {
            this.maxx = Math.max(this.maxx, cx + r);
          }
        }
        this.buf.push(`A ${r} ${r} ${0} ${0} ${side} ${x} ${y}`);
      }
      this.lastx = x;
      this.lasty = y;
    }
  };

  // src/app.ts
  var DEFAULT_CLEARANCE_PERCENT = 15;
  var clearancePercent = DEFAULT_CLEARANCE_PERCENT;
  var DEFAULT_BACKLASH_PERCENT = 0;
  var backlashPercent = DEFAULT_BACKLASH_PERCENT;
  var DEFAULT_BALANCE_PERCENT = 50;
  var balancePercent = DEFAULT_BALANCE_PERCENT;
  var DEFAUT_PRESSURE_ANGLE = 20;
  var pressureAngle = DEFAUT_PRESSURE_ANGLE;
  var DEFAULT_CONTACT_RATIO = 1.5;
  var contactRatio = DEFAULT_CONTACT_RATIO;
  var DEFAULT_PROFILE_SHIFT_PERCENT = 0;
  var profileShift = DEFAULT_PROFILE_SHIFT_PERCENT;
  var DEFAULT_GEAR_TEETH = 32;
  var gearTeeth = DEFAULT_GEAR_TEETH;
  var DEFAULT_PINION_TEETH = 16;
  var pinionTeeth = DEFAULT_PINION_TEETH;
  var DEFAULT_IS_INTERNAL = false;
  var isInternal = DEFAULT_IS_INTERNAL;
  var DEFAULT_FACE_TOL = 0.05;
  var faceTolPercent = DEFAULT_FACE_TOL;
  var DEFAULT_FILLET_TOL = 0.5;
  var filletTolPercent = DEFAULT_FILLET_TOL;
  var DEFAULT_SIZE_NUMBER = 1;
  var DEFAULT_SIZE_UNIT = "mm";
  var DEFAULT_SIZE_MEASUREMENT = "mod";
  var sizeNumber = DEFAULT_SIZE_NUMBER;
  var sizeMeasurement = DEFAULT_SIZE_MEASUREMENT;
  var sizeUnit = DEFAULT_SIZE_UNIT;
  var SIZE_UNITS = {
    px: 1,
    mm: 96 / 25.4,
    in: 96,
    ptmm: 72 / 25.4,
    ptin: 72
  };
  var animationStartTime = Date.now();
  var gearRadius = 0;
  var pinionRadius = 0;
  var stdRack = () => {
  };
  var gearRack = () => {
  };
  var pinionRack = () => {
  };
  var pinionPath;
  var gearPath;
  var pinionSvgUrl;
  var gearSvgUrl;
  var pinionSegsPerTooth = 0;
  var gearSegsPerTooth = 0;
  function update() {
    var _a;
    let params = new URLSearchParams((window.location.hash || "").substring(1));
    pressureAngle = numberParam(params, "pa", 10, 30, DEFAUT_PRESSURE_ANGLE);
    clearancePercent = numberParam(params, "clr", 0, 50, DEFAULT_CLEARANCE_PERCENT);
    backlashPercent = numberParam(params, "bkl", 0, 50, DEFAULT_BACKLASH_PERCENT);
    balancePercent = numberParam(params, "bp", 0, 100, DEFAULT_BALANCE_PERCENT);
    contactRatio = numberParam(params, "cr", 1, 2.5, DEFAULT_CONTACT_RATIO);
    profileShift = numberParam(params, "ps", -100, 100, DEFAULT_PROFILE_SHIFT_PERCENT);
    pinionTeeth = numberParam(params, "pt", 4, 200, DEFAULT_PINION_TEETH);
    gearTeeth = numberParam(params, "gt", -300, 300, DEFAULT_GEAR_TEETH);
    isInternal = false;
    if (gearTeeth < 0) {
      gearTeeth = -gearTeeth;
      isInternal = true;
    }
    if (gearTeeth < 4) {
      gearTeeth = 4;
    }
    if (isInternal && gearTeeth <= pinionTeeth) {
      gearTeeth = pinionTeeth + 1;
    }
    faceTolPercent = numberParam(params, "ft", 1e-6, 10, DEFAULT_FACE_TOL);
    filletTolPercent = numberParam(params, "fit", 1e-6, 10, DEFAULT_FILLET_TOL);
    gearRadius = gearTeeth * 0.5 / Math.PI;
    pinionRadius = pinionTeeth * 0.5 / Math.PI;
    let szLen;
    if (params.get("mod")) {
      szLen = 1 / Math.PI;
      sizeMeasurement = "mod";
    } else if (params.get("dp")) {
      szLen = 1;
      sizeMeasurement = "dp";
    } else if (params.get("cd")) {
      szLen = gearRadius + (isInternal ? -pinionRadius : pinionRadius);
      sizeMeasurement = "cd";
    } else {
      szLen = 1 / Math.PI;
      sizeMeasurement = "mod";
    }
    const szString = stringParam(params, sizeMeasurement, "");
    sizeUnit = ((_a = szString.match(/[a-zA-Z]+$/)) == null ? void 0 : _a[0]) || "";
    sizeNumber = Number(szString.substring(0, szString.length - sizeUnit.length));
    if (!isFinite(sizeNumber) || sizeNumber <= 0) {
      sizeNumber = DEFAULT_SIZE_NUMBER;
    }
    if (!SIZE_UNITS[sizeUnit]) {
      sizeUnit = "mm";
    }
    const svgScale = sizeNumber * SIZE_UNITS[sizeUnit] / szLen;
    setNumber("pa", pressureAngle);
    setNumber("cr", contactRatio);
    setNumber("gt", gearTeeth);
    setCheck("isinternal", isInternal);
    setNumber("pt", pinionTeeth);
    setNumber("ps", profileShift);
    setNumber("bp", balancePercent);
    setNumber("clr", clearancePercent);
    setNumber("bkl", backlashPercent);
    setNumber("ft", faceTolPercent);
    setNumber("fit", filletTolPercent);
    setNumber("sz", sizeNumber);
    setString("meas", sizeMeasurement);
    setString("unit", sizeUnit);
    const rackProps = {
      contactRatio,
      pressureAngle,
      profileShift,
      balancePercent,
      balanceAbsPercent: 0,
      topClrPercent: 0,
      botClrPercent: 0
    };
    stdRack = makeRack(rackProps);
    pinionRack = makeRack(__spreadProps(__spreadValues({}, rackProps), {
      botClrPercent: clearancePercent,
      balanceAbsPercent: backlashPercent * -0.5
    }));
    gearRack = makeRack(__spreadProps(__spreadValues({}, rackProps), {
      balancePercent: isInternal ? rackProps.balancePercent : 100 - rackProps.balancePercent,
      profileShift: isInternal ? profileShift : -profileShift,
      botClrPercent: isInternal ? 0 : clearancePercent,
      topClrPercent: isInternal ? clearancePercent : 0,
      balanceAbsPercent: backlashPercent * (isInternal ? 0.5 : -0.5)
    }));
    const faceT = faceTolPercent / (100 * Math.PI);
    const filletT = filletTolPercent / (100 * Math.PI);
    const pinionCutter = new GearCutter(pinionTeeth, pinionRadius, faceT, filletT);
    pinionRack(new XFormPen(pinionCutter).rotate(-90).translate(0, pinionRadius), true);
    const pinionRecorder = new RecordingPen();
    pinionCutter.drawToothPath(pinionRecorder, true);
    pinionPath = pinionRecorder.path;
    pinionSegsPerTooth = pinionRecorder.countSegments();
    const gearCutter = new GearCutter(gearTeeth, gearRadius, faceT, filletT);
    gearRack(new XFormPen(gearCutter).rotate(-90).translate(0, gearRadius), true);
    const gearRecorder = new RecordingPen();
    gearCutter.drawToothPath(gearRecorder, true);
    gearPath = gearRecorder.path;
    gearSegsPerTooth = gearRecorder.countSegments();
    if (pinionSvgUrl) {
      URL.revokeObjectURL(pinionSvgUrl);
    }
    pinionSvgUrl = createSvgObjectUrl(svgScale, pinionPath, pinionRadius, pinionTeeth);
    document.getElementById("img1").src = pinionSvgUrl;
    document.getElementById("link1").href = pinionSvgUrl;
    document.getElementById("link1").download = `pinion${pinionTeeth}pa${pressureAngle}.svg`;
    document.getElementById("svgTitle1").innerText = `Pinion (${pinionSegsPerTooth} arcs/tooth)`;
    if (gearSvgUrl) {
      URL.revokeObjectURL(gearSvgUrl);
    }
    gearSvgUrl = createSvgObjectUrl(svgScale, gearPath, gearRadius, gearTeeth);
    document.getElementById("img2").src = gearSvgUrl;
    document.getElementById("link2").href = gearSvgUrl;
    document.getElementById("link2").download = `gear${gearTeeth}pa${pressureAngle}.svg`;
    document.getElementById("svgTitle2").innerText = `Gear (${gearSegsPerTooth} arcs/tooth)`;
  }
  function submit() {
    const fields = [
      ["pa", DEFAUT_PRESSURE_ANGLE],
      ["cr", DEFAULT_CONTACT_RATIO],
      ["gt", -1],
      ["pt", -1],
      ["ps", DEFAULT_PROFILE_SHIFT_PERCENT],
      ["bp", DEFAULT_BALANCE_PERCENT],
      ["clr", DEFAULT_CLEARANCE_PERCENT],
      ["bkl", DEFAULT_BACKLASH_PERCENT],
      ["ft", DEFAULT_FACE_TOL],
      ["fit", DEFAULT_FILLET_TOL]
    ];
    let parts = [];
    for (const [id, defval] of fields) {
      let val = document.getElementById(id).value.trim();
      if (id === "gt" && document.getElementById("isinternal").checked) {
        val = "-" + val;
      }
      if (val !== String(defval)) {
        parts.push(id + "=" + val);
      }
    }
    const sz = Number(document.getElementById("sz").value) || DEFAULT_SIZE_NUMBER;
    const meas = document.getElementById("meas").value;
    const unit = document.getElementById("unit").value;
    if (meas !== DEFAULT_SIZE_MEASUREMENT || sz !== DEFAULT_SIZE_NUMBER || unit !== DEFAULT_SIZE_UNIT) {
      parts.push(`${meas}=${sz}${unit}`);
    }
    let newhash = "#" + parts.join("&");
    if (newhash != window.location.hash) {
      window.location.hash = newhash;
    } else {
      update();
    }
  }
  function createSvgObjectUrl(svgScale, path, pitchRadius, nTeeth) {
    const svg = new SvgRecorder({
      scale: -svgScale
    });
    svg.drawCircle(
      {
        stroke: "blue",
        strokeWidth: 0.01
      },
      0,
      0,
      pitchRadius
    );
    svg.draw(
      {
        stroke: "black",
        strokeWidth: 0.01
      },
      (pen, domove) => {
        for (let i = 0; i < nTeeth; ++i) {
          const p = new XFormPen(pen);
          p.rotate(i * 360 / nTeeth);
          path(p, domove);
          domove = false;
        }
      }
    );
    const blob = new Blob([svg.getSvgText()], { type: "image/svg+xml" });
    return URL.createObjectURL(blob);
  }
  function setNumber(id, val) {
    document.getElementById(id).value = String(val);
  }
  function setString(id, val) {
    document.getElementById(id).value = val;
  }
  function setCheck(id, val) {
    document.getElementById(id).checked = val;
  }
  function numberParam(params, name, minval, maxval, defval) {
    let v = params.get(name);
    if (v == null) {
      return defval;
    }
    v = Number(v);
    if (isNaN(v)) {
      return defval;
    }
    if (minval != null && v < minval) {
      return minval;
    }
    if (maxval != null && v > maxval) {
      return maxval;
    }
    return v;
  }
  function stringParam(params, name, defval) {
    return params.get(name) || defval;
  }
  function drawRack(pen, rack, shift, mint, maxt) {
    mint -= shift;
    maxt -= shift;
    let start = Math.floor(mint + 0.5);
    const xpen = new XFormPen(pen);
    rack(xpen.copy().translate(start + shift, 0), true);
    for (let t = start + 1; t < maxt + 0.9; t++) {
      rack(xpen.copy().translate(t + shift, 0), false);
    }
  }
  function animationFrame() {
    let canvas = document.getElementById("canvas");
    const cwid = canvas.width;
    const chei = canvas.height;
    let ctx = canvas.getContext("2d");
    ctx.save();
    ctx.clearRect(0, 0, cwid, chei);
    let scale = cwid / 3;
    const pen = new XFormPen(new CanvasPen(ctx)).translate(cwid / 2, chei / 2).scale(scale, true);
    let shift = Date.now() - animationStartTime;
    shift /= 4e3;
    shift -= Math.floor(shift);
    ctx.beginPath();
    ctx.strokeStyle = "#A0A0A0";
    drawRack(pen, stdRack, shift, -1.5, 1.5);
    ctx.stroke();
    ctx.strokeStyle = "#000000";
    const bkwAdjust = backlashPercent / (400 * Math.PI);
    if (pinionPath) {
      ctx.beginPath();
      const p = pen.copy().translate(0, -pinionRadius).rotate(90 - 360 / pinionTeeth);
      p.rotate((shift + bkwAdjust) * -360 / pinionTeeth);
      pinionPath(p, true);
      p.rotate(360 / pinionTeeth);
      pinionPath(p, true);
      p.rotate(360 / pinionTeeth);
      pinionPath(p, true);
      p.rotate(360 / pinionTeeth);
      pinionPath(p, true);
      p.rotate(360 / pinionTeeth);
      pinionPath(p, true);
      ctx.stroke();
    }
    if (gearPath) {
      ctx.beginPath();
      let p;
      let sh = shift;
      if (isInternal) {
        p = pen.copy().translate(0, -gearRadius).rotate(90 - 360 / gearTeeth);
      } else {
        p = pen.copy().translate(0, gearRadius).scale(1, true).rotate(90 - 360 / gearTeeth);
        sh += -0.5;
        sh -= Math.floor(sh);
      }
      p.rotate((sh - bkwAdjust) * -360 / gearTeeth);
      gearPath(p, true);
      p.rotate(360 / gearTeeth);
      gearPath(p, true);
      p.rotate(360 / gearTeeth);
      gearPath(p, true);
      p.rotate(360 / gearTeeth);
      gearPath(p, true);
      p.rotate(360 / gearTeeth);
      gearPath(p, true);
      ctx.stroke();
    }
    ctx.restore();
    requestAnimationFrame(animationFrame);
  }
  update();
  window.onhashchange = update;
  document.getElementById("submit").onclick = submit;
  requestAnimationFrame(animationFrame);
})();
