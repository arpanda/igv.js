import t_dist from './t_dist.js'
import baseCNVpytorVCF from './baseCNVpytorVCF.js'


// TODO -- remove this hardcoded value
const genome_size = 2871000000;

class MeanShiftCaller extends baseCNVpytorVCF{
    /**
     * Creates an instance of CombinedCaller.
     * 
     * @param {Array} wigFeatures - An array of arrays containing wig formatted data for each chromosome and bin.
     * @param {number} binSize - The size of the bins used in the wig data.
     * @param {string} refGenome - reference genome name
     */
    constructor(wigFeatures, binSize, refGenome) {
        super(wigFeatures, binSize, refGenome)
        this.binBands = [2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 112, 128]
        
    }

    async callMeanshift(repeats = 3){
        // applying gc correction
        await this.apply_gcCorrection()
        console.log('wigFeatures', this.wigFeatures)
        
        let partitionLevels = await this.partition(repeats=repeats)
        console.log("Partition: ", partitionLevels)

        Object.entries(this.wigFeatures).forEach(([chr, chrRD]) => {
            chrRD.forEach((bin, index) => {
                if (partitionLevels[chr]){ 
                    bin.partitionLevel = parseInt(partitionLevels[chr][index])
                }
            });
        })
        // console.log("WigFeatures: ", this.wigFeatures)

        let cnvs = await this.cnvCalling(partitionLevels)
        console.log('cnvs: ', cnvs)
        
        this.CNVcalls = {};
        this.CNVcalls[this.binSize] = {}
        this.CNVcalls[this.binSize]['ReadDepth'] = cnvs


        var rawbinScore = this.formatDataStructure('binScore', this.globalMean)
        var gcCorrectedBinScore = this.formatDataStructure('gcCorrectedBinScore', this.globalMean)
        var partitionBinScore = this.formatDataStructure('partitionLevel', this.globalMean)
       
        const fetchedData = {binScore: rawbinScore, gcCorrectedBinScore: gcCorrectedBinScore, segmentsCNV: partitionBinScore}
        return fetchedData
        

    }
    
    getRDSignalBandWidth(data_array) {
        const threshold = this.globalMean / 4;

        // The values are reversed and squared as they will be used to calculate a gradient function.
        // This optimization is done to speed up the calculations.

        const constantValue = 4 / this.globalStd ** 2;
        return data_array.map(value => {
            return value > threshold ? this.globalMean / (this.globalStd ** 2 * value) : constantValue;
        });
    }

    async partition(repeats = 3){
        
        // sort the dictionary based on chromosome names;
        let sortedDictionary = {};
        Object.keys(this.wigFeatures).sort((a, b) => a.localeCompare(b, undefined, {numeric: true})).forEach(key => {
            sortedDictionary[key] = this.wigFeatures[key];
        });

        let binScoreField = this.gcFlag ? "gcCorrectedBinScore": "binScore" ;

        var chrLevels = {}
        // Object.entries(this.wigFeatures).forEach(([chr, chr_rd]) => {
        for (const [chr, chrWig] of Object.entries(sortedDictionary)) {

            // boolean array; Initiate with all false values
            var masked = new Array(chrWig.length).fill(false)

            // set the level; score from either RAW or GC corrected bin score 
            var levels = chrWig.map((item, index) => !masked[index] ? item[binScoreField] : undefined);
            
            this.binBands.forEach((bin_band, bin_band_index) => {
                
                // not masked levels at current bin
                // get boolean values
                var not_masked = masked.map((value, index) => { return !value; })

                // not mask level
                var nm_levels = levels.filter((_, index) => !masked[index]);
                
                // set the mask border
                var mask_borders = [0]
                var count = 0
                
                // the masked array was declared previously
                masked.forEach(item => {
                    if (item) {
                        if (count > 0) {
                            mask_borders.push(mask_borders[mask_borders.length - 1] + count - 1);
                            count = 0;
                        }
                    } else if (!item) { count++; }
                });

                mask_borders.shift()

                // repeating steps
                for (let step = 0; step < repeats; step++) {
                    var isig = this.getRDSignalBandWidth(nm_levels)
                    
                    // get the direction of meanshift vector for a bin
                    // it compares near by bins to get the direction of the vector; bin band defines the range of comparison

                    var grad = new Array(nm_levels.length).fill(0);

                    for (let i = 0; i < nm_levels.length; i++) {
                        const start_bin = Math.max(0, i - 3 * bin_band);
                        const end_bin = Math.min(nm_levels.length - 1, i + 3 * bin_band + 1); 
                        for (let j = start_bin ; j <= end_bin; j++) {

                            var g_value = (j - i) * Math.exp((-0.5 * (j - i) ** 2) / bin_band ** 2) *  Math.exp(-0.5 * (nm_levels[i] - nm_levels[j]) ** 2 * isig[i]);

                            grad[i] += g_value
                        }
                    }
                    
                    // get the border; if there is a change of gradient, it is a border
                    var border = new Array();
                    for (var i = 0; i < grad.length - 1; i++) {
                        if ((grad[i] < 0) & (grad[i + 1] >= 0)) border.push(i);
                    }

                    border.push(grad.length - 1)
                    border = border.concat(mask_borders).sort((a, b) => a - b)
                    border = Array.from(new Set(border))

                    var pb = 0;
                    for (var i = 0; i < border.length; i++) {
                        var range_array = nm_levels.slice(pb, border[i] + 1)
                        var range_mean = range_array.reduce((acc, n) => acc + n) / range_array.length

                        nm_levels.fill(range_mean, pb, border[i] + 1)
                        pb = border[i] + 1
                    }
                }

                for (var i = 0, j = 0; i < levels.length; i++) {
                    if (not_masked[i]) {
                        levels[i] = nm_levels[j]
                        j++
                    }
                }

                //get the border
                var border = new Array();
                for (var i = 0; i < levels.length - 1; i++) {
                    //if(i == levels.length -1) continue;
                    var diff = Math.abs(levels[i + 1] - levels[i]);

                    if (diff > 0.01) border.push(i + 1);
                }

                border.unshift(0);
                border.push(levels.length);

                // reset the mask
                masked = new Array(this.wigFeatures.length).fill(false);

                // check the borders
                for (var i = 1; i < border.length; i++) {
                    var seg = [border[i - 1], border[i]]
                    var seg_left = [border[i - 1], border[i - 1]]
                    if (i > 1) { seg_left[0] = border[i - 2] } else continue;

                    var seg_right = [border[i], border[i]];
                    if (i < border.length - 1) { seg_right[1] = border[i + 1] } else continue;

                    var n = seg[1] - seg[0];
                    var n_left = seg_left[1] - seg_left[0];
                    var n_right = seg_right[1] - seg_right[0];
                    if (n <= 1) continue;
                    var seg_array = new DataStat(levels.slice(seg[0], seg[1]));

                    if (n_right <= 15 || n_left <= 15 || n <= 15) {
                        var ns = 1.8 * Math.sqrt(levels[seg_left[0]] / this.globalMean) * this.globalStd;
                        if (Math.abs(levels[seg_left[0]] - levels[seg[0]]) < ns) { continue }

                        ns = 1.8 * Math.sqrt(levels[seg_right[0]] / this.globalMean) * this.globalStd;
                        if (Math.abs(levels[seg_right[0]] - levels[seg[0]]) < ns) { continue }
                    } else {
                        var seg_left_array = levels.slice(seg_left[0], seg_left[1])
                        var seg_left_1 = new DataStat(seg_left_array);

                        var seg_right_array = levels.slice(seg_right[0], seg_right[1])
                        var seg_right_1 = new DataStat(seg_right_array);

                        var ttest_2sample_1 = t_test_2_samples(seg_array.mean, seg_array.std, seg_array.data.length,
                            seg_left_1.mean, seg_left_1.std, seg_left_1.data.length);
                        if (ttest_2sample_1 > (0.01 / genome_size) * this.binSize * (n + n_left)) { continue }

                        var ttest_2sample_2 = t_test_2_samples(seg_array.mean, seg_array.std, seg_array.data.length,
                            seg_right_1.mean, seg_right_1.std, seg_right_1.data.length);
                        if (ttest_2sample_2 > (0.01 / genome_size) * this.binSize * (n + n_right)) { continue }
                    }

                    var ttest_1sample_1 = t_test_1_sample(this.globalMean, seg_array.mean, seg_array.std, seg_array.data.length)
                    if (ttest_1sample_1 > 0.05) { continue }
                    let segments_rd = nm_levels.slice(seg[0], seg[1])
                   
                    var raw_seg_data = new DataStat(segments_rd);

                    masked.fill(true, seg[0], seg[1]);
                    levels.fill(raw_seg_data.mean, seg[0], seg[1]);
                }
        
            });
            // console.log("after applying partition: ", levels)
            chrLevels[chr] = levels
            // break
        }
        
        return chrLevels
    }

    async cnvCalling(levels) {

        var delta = 0.25 * this.globalMean
        var min = this.globalMean - delta, max = this.globalMean + delta;
        var normal_genome_size = 2971000000

        // var levels = this.meanShiftCaller(bin_size)
        // var levels = this.MeanShiftCallerV2(bin_size)
        
        // var cnv_levels = {}   
        var cnvCalls = {}
        Object.entries(levels).forEach(([chr, chr_levels]) => {
            
            var done = false
            while (!done) {
                done = true

                //
                // get the borders
                //
                var borders = new Array(1).fill(0);
                for (let i = 0; i < chr_levels.length - 1; i++) {
                    var diff = Math.abs(chr_levels[i + 1] - chr_levels[i]);
                    if (diff > 0.01) borders.push(i + 1);
                }
                // console.log('chr_levels', chr_levels)
                borders.push(chr_levels.length);

                for (let ix = 0; ix < borders.length - 2; ix++) {
                    var v1 = Math.abs(chr_levels[borders[ix]] - chr_levels[borders[ix + 1]])
                    
                    if (v1 < delta) {
                        var v2 = v1 + 1, v3 = v1 + 1;

                        if (ix > 0) { v2 = Math.abs(chr_levels[borders[ix]] - chr_levels[borders[ix - 1]]) }
                        if (ix < borders.length - 3) { v3 = Math.abs(levels[borders[ix + 1]] - chr_levels[borders[ix + 2]]) }

                        if (v1 < v2 && v1 < v3) {
                            done = false

                            var tmp_array = new DataStat(chr_levels.slice(borders[ix], borders[ix + 2]))
                            chr_levels.fill(tmp_array.mean, borders[ix], borders[ix + 2]);
                            borders.splice(ix + 1, ix + 1);
                        }
                    }
                }
            }

            // var chr_rd = this.rd[chr]
            var chr_rd = []
            // Object.entries(this.wigFeatures[chr]).forEach(([bin, binDict]) => { chr_rd.push(binDict.binScore) });
            Object.entries(this.wigFeatures[chr]).forEach(([bin, binDict]) => { chr_rd.push(binDict.partitionLevel) });
            // console.log(chr, chr_rd)

            //
            // Calling Segments
            //

            //var flags = [""] * levels.length;
            var flags = new Array(chr_levels.length).fill("")
            var segments = []

            var b = 0
            var pval = (0.05 * this.binSize) / normal_genome_size
            while (b < chr_levels.length) {
                var b0 = b, border_start = b;
                while ((b < chr_levels.length) & (chr_levels[b] < min)) {
                    b += 1
                }
                var border_end = b
               
                if (border_end > border_start + 1) {
                    // console.log('border min', border_start, border_end)
                    var adj = adjustToEvalue(this.globalMean, this.std, chr_rd, border_start, border_end, pval)
                    if (adj) {
                        border_start, border_end = adj;
                        segments.push([border_start, border_end + 1])
                        flags.fill("D", border_start, border_end)
                    }
                }

                border_start = b;
                while ((b < chr_levels.length) & (chr_levels[b] > max)) { b += 1 }
                border_end = b
                
                if (border_end > border_start + 1) {
                    adj = adjustToEvalue(this.globalMean, this.std, chr_rd, border_start, border_end, pval)
                    
                    if (adj) {
                        border_start, (border_end = adj);
                        segments.push([border_start, border_end, +1]);
                        // flags[bs:be] = ["A"] * (be - bs)
                        flags.fill("A", border_start, border_end);
                    }
                }
                if (b == b0) b += 1;
            }

            //
            //  Calling additional deletions
            //
            b = 0;
            while (b < chr_levels.length) {
                while ((b < chr_levels.length) & (flags[b] != "")) b += 1;
                border_start = b;
                while ((b < chr_levels.length) & (chr_levels[b] < min)) b += 1;
                border_end = b;
                if (border_end > border_start + 1) {
                    if (gaussianEValue(this.globalMean, this.std, chr_rd, border_start, border_end) < 0.05 / normal_genome_size) {
                        segments.push([border_start, border_end, -1])
                        flags.fill(["d"] * (border_end - border_start), border_start, border_end)
                    }
                    b -= 1;
                }
                b += 1;
            }

            b = 0;
            var cf;
            if (b < chr_levels.length) cf = flags[b];

            border_start = 0;

            //var merge = [...this.rd]
            var merge = [...chr_rd]
            // console.log('initial merge', merge)
            while (b < chr_levels.length) {
                while (flags[b] == cf) {
                    b += 1;
                    if (b >= flags.length) break;
                }
                if (b > border_start) {
                    var merge_arr = new DataStat(merge.slice(border_start, b));
                    merge.fill(merge_arr.mean, border_start, b);
                }
                if (b < chr_levels.length) cf = flags[b];
                border_start = b;
            }

            // calculate calls
            b = 0
            while (b < chr_levels.length) {
                cf = flags[b]
                if (cf == "") {
                    b += 1
                    continue
                }
                border_start = b
                while (b < chr_levels.length & cf == flags[b]) { b += 1 }
                var border_arr = new DataStat(merge.slice(border_start, b))
                let normalized_cn = border_arr.mean / this.globalMean

                let cnvType = cf === "D" ? "deletion" : "duplication";
                
                let cnv_dict = {
                    chr: chr,
                    start: this.binSize * border_start + 1,
                    end: this.binSize * b,
                    size: this.binSize * (b - border_start + 1),
                    cn: normalized_cn ,
                    type: cnvType,
                    bins: (b - border_start + 1)
                }
                if (!cnvCalls[cnvType]) { cnvCalls[cnvType] = [] }
                cnvCalls[cnvType].push(cnv_dict)

            }
        });

        // return cnv_levels
        return cnvCalls
    }

}

class DataStat {
    constructor(data_array) {
        this.data = data_array
        this.mean = data_array.reduce((acc, n) => acc + n) / data_array.length
        this.std = Math.sqrt(data_array.reduce((acc, n) => (n - this.mean) ** 2) / data_array.length)
    }
}

function erf(x) {
    var m = 1.0, s = 1.0, sum = x * 1.0;
    for (var i = 1; i < 50; i++) {
        m *= i;
        s *= -1;
        sum += (s * Math.pow(x, 2.0 * i + 1.0)) / (m * (2.0 * i + 1.0));
    }
    return (2 * sum) / Math.sqrt(3.14159265358979)
}

function getEValue(mean, sigma, rd, start, end) {
    var arr = new DataStat(rd.slice(start, end));
    
    if (arr.std == 0) {
        if (sigma > 0) { arr.std = (sigma * arr.mean) / mean }
        else { arr.std = 1 }
    }
    var p_val = t_test_1_sample(mean, arr.mean, arr.std, end - start) / (end - start)
    return p_val
}

function gaussianEValue(mean, sigma, rd, start, end) {
    var arr = new DataStat(rd.slice(start, end))

    if (arr.mean < mean) {
        var x = (arr.max - arr.mean) / (sigma * Math.sqrt(2))
        return Math.pow(0.5 * (1 + erf(x)), end - start)
    }
    var x = (arr.min - arr.mean) / (sigma * Math.sqrt(2))
    return Math.pow(0.5 * (1 - erf(x)), end - start)
}

function adjustToEvalue(mean, sigma, rd, start, end, pval, max_steps = 1000) {
    var val = getEValue(mean, sigma, rd, start, end)

    var step = 0, done = false
    while ((val > pval) & !done & (step < max_steps)) {
        done = true
        step += 1
        var [v1, v2, v3, v4] = [1e10, 1e10, 1e10, 1e10];
        if (start > 0) v1 = getEValue(mean, sigma, rd, start - 1, end);
        if (end - start > 2) {
            var v2 = getEValue(mean, sigma, rd, start + 1, end)
            var v3 = getEValue(mean, sigma, rd, start, end - 1)
        }
        if (end < rd.length) { var v4 = getEValue(mean, sigma, rd, start, end + 1) }
        if (Math.min[(v1, v2, v3, v4)] < val) {
            done = false
            if (v1 == Math.min[(v1, v2, v3, v4)]) {
                start -= 1;
                val = v1;
            }
            elif(v2 == Math.min[(v1, v2, v3, v4)]); {
                start += 1;
                val = v2;
            }
            elif(v3 == Math.min[(v1, v2, v3, v4)]); {
                end -= 1;
                val = v3;
            }
            elif(v4 == Math.min[(v1, v2, v3, v4)]); {
                end += 1;
                val = v4;
            }
        }
    }
    if (val <= pval) { return start, end }
    return 0;
}

function t_test_1_sample(mean, m, s, n) {
    if (s == 0) s = 1;
    var t = ((mean - m) / s) * Math.sqrt(n)
    var p = 1.0 - t_dist.TdistributionCDF(Math.abs(t), (n - 1))
    return p
}

function t_test_2_samples(m1, s1, n1, m2, s2, n2) {
    if (s1 == 0) s1 = 1;
    if (s2 == 0) s2 = 1;
    var t = (m1 - m2) / Math.sqrt(s1 ** 2 / n1 + s2 ** 2 / n2);
    var df = ((s1 ** 2 / n1 + s2 ** 2 / n2) ** 2 * (n1 - 1) * (n2 - 1)) /
        ((s1 ** 4 * (n2 - 1)) / n1 ** 2 + (s2 ** 4 * (n1 - 1)) / n2 ** 2);

    var p = 1.0 - t_dist.TdistributionCDF(Math.abs(t), parseInt(df + 0.5))

    return p
}


export default { MeanShiftCaller };
