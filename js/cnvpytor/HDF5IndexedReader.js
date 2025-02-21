import {openH5File} from "../../node_modules/hdf5-indexed-reader/dist/hdf5-indexed-reader.esm.js"


class SignalNames{
    /**
     * 
     * @param {string} chrom - chromosome name
     * @param {integer} bin_size - bin size
     */
    constructor(chrom, bin_size){
        this.chrom = chrom
        this.signal_bin_size = bin_size

        let rd_flag = ""
        this.signals = {
            'RD_Raw': `his_rd_p_${this.chrom}_${this.signal_bin_size}${rd_flag}`,
            'RD_Raw_gc_corrected': `his_rd_p_${this.chrom}_${this.signal_bin_size}_GC`,
            'rd_partition' : `his_rd_p_${this.chrom}_${this.signal_bin_size}_partition_GC_merge`,
            'baf': `snp_likelihood_${this.chrom}_${this.signal_bin_size}_mask`,
            'baf_i1': `snp_i1_${this.chrom}_${this.signal_bin_size}_mask`,
            'Mosaic_segments' : `his_rd_p_${this.chrom}_${this.signal_bin_size}_partition_GC_mosaic_segments_2d`,
            'Mosaic_calls': `his_rd_p_${this.chrom}_${this.signal_bin_size}_partition_GC_mosaic_call_2d`,
            'ReadDepth_calls': `calls_${this.chrom}_${this.signal_bin_size}_GC`,
            "2D_calls": `calls_2d_${this.chrom}_${this.signal_bin_size}_mask_GC_mask`
        }
    }
}


class HDF5Reader {
    /**
     * 
     * @param {string} h5_file - path for the pytor file
     * @param {integer} bin_size - bin size
     */
    constructor(h5_file, bin_size=100000){

        this.h5_file = h5_file;
        this.bin_size = bin_size;
        this.h5_obj = undefined
        this.pytorKeys = [];
        this.availableBins = [];
        this.processedChroms = {};
        this.allWigFeatures = {};
        this.CNVcalls = {};
        this.rd_stats = {}
    }
    
    async fetch(){

        if(!this.h5_obj) {
            this.h5_obj = await openH5File({
                url: this.h5_file,
                fetchSize: 1000000,
                maxSize: 200000000
            })
        }
        return this.h5_obj
    }

    /**
     * 
     * @returns - a list of keys of the pytor file
     */
    async get_keys(){
        let h5_obj = await this.fetch();
        return h5_obj.keys
    }

    async get_rd_signal(bin_size = this.bin_size, cnv_caller="2D", chrom=undefined){
        // Fetch the pytor file and get keys
        const h5Obj = await this.fetch();
        this.pytorKeys = h5Obj.keys;

        // get available bin sizes
        const signalBin = new ParseSignals(this.pytorKeys);
        this.availableBins = signalBin.getAllBins();
        
        // check if the user provided bin is available, else set the last bin_size
        if(! this.availableBins.includes(bin_size)){
            bin_size = this.availableBins.at(-1);    
        }
        
        // get rd chromosomes and rd stat
        const rdChromosomes = await this.getChromosomes(chrom);
        if (!this.rd_stats[bin_size]){
            this.rd_stats[bin_size] = await this.get_rd_stat(bin_size)
        }

        // prepare wig formatted file for all chromosome
        
        // save the chroms and bin_sizes that are already processed
        if (!this.allWigFeatures[bin_size]){
            this.processedChroms[bin_size] = rdChromosomes
        }else{
            this.processedChroms[bin_size].push(rdChromosomes)
        }
        
        // get signals; two types of signal, 1: common 2: callers
        await this.getCommonWigfeatures(rdChromosomes, bin_size)
        await this.getCallerWigfeatures(cnv_caller, rdChromosomes, bin_size)

        // set callers name
        this.getCallersName()
        
    }
    getCallersName(){
        const patterns = [
            { name: "ReadDepth", pattern: /calls_[\w\d]+_(\d+)_GC/ },
            { name: "2D", pattern: /calls_2d_[\w\d]+_(\d+)_mask_GC_mask/ },
        ]
        let callers = {}
        for (const { name, pattern } of patterns) {
            callers[name] = []; // Initialize an empty array for each pattern
    
            this.pytorKeys.forEach(item => {
                const match = item.match(pattern);
                if (match) {
                    callers[name].push(parseInt(match[1])); // Collect all matching bin sizes
                }
            });
    
            // Remove duplicates if necessary
            callers[name] = [...new Set(callers[name])]; 
        }
        this.callers = Object.keys(callers);
        return callers
    }


    async getCommonWigfeatures(rdChromosomes, binSize){
        
        if (!this.allWigFeatures[binSize]){
            this.allWigFeatures[binSize] = {}
        }
        let rdStat = this.rd_stats[binSize]
        let rd_signals = ['RD_Raw',  'RD_Raw_gc_corrected']

        for (const chrom of rdChromosomes) {
            const signalNameObj = new SignalNames(chrom, binSize);
            
            // processing rd signals 
            for (const rd_signal of rd_signals){
                let signal_name = signalNameObj.signals[rd_signal]
                if (this.pytorKeys.includes(signal_name)){
                    let chr_wig = await this.get_chr_signal(chrom, binSize, signal_name, rdStat)
                    if(!this.allWigFeatures[binSize][rd_signal]){
                        this.allWigFeatures[binSize][rd_signal] = []
                    }
                    this.allWigFeatures[binSize][rd_signal].push(...chr_wig)
                }
            }
            
            const [baf1, baf2] = await this.getBafSignals(chrom, binSize, signalNameObj.signals.baf_i1);
            let baf_signal = 'BAF1'
            if (!this.allWigFeatures[binSize][baf_signal]){
                this.allWigFeatures[binSize][baf_signal] = []
            }
            this.allWigFeatures[binSize][baf_signal].push(...baf1);

            baf_signal = 'BAF2'
            if (!this.allWigFeatures[binSize][baf_signal]){
                this.allWigFeatures[binSize][baf_signal] = []
            }
            this.allWigFeatures[binSize][baf_signal].push(...baf2);
            
        }
        
    }
    async getCallerWigfeatures(cnv_caller, rdChromosomes, binSize){
        let rdStat = this.rd_stats[binSize]
        
        if (!this.CNVcalls[binSize]){
            this.CNVcalls[binSize] = {}
        }
        if (!this.CNVcalls[binSize][cnv_caller]){
            this.CNVcalls[binSize][cnv_caller] = {}
        }
        for (const chrom of rdChromosomes) {
            const signalNameObj = new SignalNames(chrom, binSize);

            if (!this.allWigFeatures[binSize][cnv_caller]){
                this.allWigFeatures[binSize][cnv_caller] = []
            }
    
            if (cnv_caller == '2D'){   
                this.allWigFeatures[binSize][cnv_caller].push(...await this.rd_call_combined(chrom, binSize, rdStat, signalNameObj));

            }else if(cnv_caller == 'ReadDepth'){
                this.allWigFeatures[binSize][cnv_caller].push(...await this.get_chr_signal(chrom, binSize, signalNameObj.signals.rd_partition, rdStat));
                
            }
            
            let signal_name = `${cnv_caller}_calls`
            // console.log(signalNameObj.signals)
            
            let calls = await this.get_calls(chrom, signalNameObj.signals[signal_name], cnv_caller)
            // console.log(calls)
            if(calls){
                for (const [cnv_type, cnv_list] of Object.entries(calls)) {
                    if(!this.CNVcalls[binSize][cnv_caller][cnv_type]){
                        this.CNVcalls[binSize][cnv_caller][cnv_type] = cnv_list
                    }else{
                        this.CNVcalls[binSize][cnv_caller][cnv_type].push(...cnv_list)
                    }
                }
            }

        }
    }

    async getWigFeatures(rdChromosomes, binSize, rdStat) {
        const wigFeatures = {
            RD_Raw: [],
            RD_Raw_gc_coor: [],
            ReadDepth: [],
            "2D": [],
            ReadDepthCall: {},
            BAF1: [],
            BAF2: []
        };

        for (const chrom of rdChromosomes) {
            const signalNameObj = new SignalNames(chrom, binSize);

            wigFeatures.RD_Raw.push(...await this.get_chr_signal(chrom, binSize, signalNameObj.signals.RD_Raw, rdStat));
            wigFeatures.RD_Raw_gc_coor.push(...await this.get_chr_signal(chrom, binSize, signalNameObj.signals.RD_Raw_gc_corrected, rdStat));
            wigFeatures.ReadDepth.push(...await this.get_chr_signal(chrom, binSize, signalNameObj.signals.rd_partition, rdStat));

            wigFeatures["2D"].push(...await this.rd_call_combined(chrom, binSize, rdStat, signalNameObj));

            const [baf1, baf2] = await this.getBafSignals(chrom, binSize, signalNameObj.signals.baf_i1);

            // 2d calls
            let calls = await this.get_calls(chrom, signalNameObj.signals['2Dcalls'], "2D")
            // console.log("calls: ", calls)
            
            // ReadDepth_calls
            // let calls = await this.get_calls(chrom, signalNameObj.signals.ReadDepth_calls, "ReadDepth")
            // console.log(calls)
            if(calls){
                
                for (const [cnv_type, cnv_list] of Object.entries(calls)) {
                    if(!wigFeatures['ReadDepthCall'][cnv_type]){
                        wigFeatures['ReadDepthCall'][cnv_type] = cnv_list
                    }else{
                        wigFeatures['ReadDepthCall'][cnv_type].push(...cnv_list)
                    }
                }
                
            }

            wigFeatures.BAF1.push(...baf1);
            wigFeatures.BAF2.push(...baf2);
        }
        return wigFeatures;
    }

    async getChromosomes(refChroms) {
        // return chromosome names if they exists in the rd_chromosomes
        const rdChroms_obj = await this.h5_obj.get("rd_chromosomes");
        const rdChroms = await rdChroms_obj.value
        if(!refChroms){
            return rdChroms 
        }else{
            let refChromsSet = new Set(refChroms)
            return rdChroms.filter(item => refChromsSet.has(item));
            
        }
    }

    setCallers(wigFeatures) {
        this.callers = [];
        if (wigFeatures.ReadDepth.length) this.callers.push('ReadDepth');
        if (wigFeatures["2D"].length) this.callers.push('2D');
    }

    decode_segments(segments_arr){
        let max = 2 ** 32 - 1
        let segments = []
        let l = []
        for (let x of segments_arr){
            if(x == max){
                segments.push(l)
                l = []
            } else{
                l.push(x)
            }
        }
        return segments
    }

    async get_calls(chrom, signal_name, caller="2D"){

        if (this.pytorKeys.includes(signal_name)){
            const chrom_dataset = await this.h5_obj.get(signal_name)
            let chrom_data = await chrom_dataset.to_array();
            
            let typeMap = new Map([ [-1, "deletion"], [0, "cnnloh"], [1, "duplication"]]);
            
            let keys_meanshift = ["start", "end", "size", "cnv", "p_val", "p_val_2", "p_val_3", "p_val_4", "Q0", "pN", "dG"];
    
            let keys_2d = ["start", "end", "size", "cnv", "p_val", "lh_del", "lh_loh", "lh_dup", "Q0", "bins", "baf", "rd_p_val", "baf_p_val", "segment", "hets", "homs", "pN", "pNS", "pP"]
            let keys = caller === '2D'? keys_2d:keys_meanshift
            let calls = {}
            for (let call of chrom_data){
                let cnvType = typeMap.get(call[1])
                let row_info = Object.fromEntries(call.slice(2).map((_, idx) => [keys[idx], _]));
                row_info['chr'] = chrom
                if (!calls[cnvType]) { calls[cnvType] = [] }
                calls[cnvType].push(row_info)
            }

            return calls
        }
    }
    async rd_call_combined(chrom, bin_size, rd_stat, signal_name_obj){
        let chr_wig = [];
        
        let segments
        let mosaic_call_segments = signal_name_obj.signals['Mosaic_segments']
        if (this.pytorKeys.includes(mosaic_call_segments)){
            const chrom_dataset = await this.h5_obj.get(mosaic_call_segments)
            let chrom_data = await chrom_dataset.value
            segments = this.decode_segments(chrom_data)
            
        }

        let mosaic_calls = signal_name_obj.signals['Mosaic_calls']
        if (this.pytorKeys.includes(mosaic_calls)){
            const segments_call_dataset = await this.h5_obj.get(mosaic_calls)
            let segments_call = await segments_call_dataset.to_array() //create_nested_array(value, shape)
            segments.forEach((ind_segment, segment_idx) => {
                ind_segment.forEach((bin_value, bin_idx) => {
                    chr_wig.push({chr:chrom, start: bin_value*bin_size, end: (bin_value+1) * bin_size, value: (segments_call[0][segment_idx]/rd_stat[4]) *2})
                })
            })
        }

        return chr_wig
        
    }
    
    /**
     * returns a list for rd statistics information 
     * @param {integer} bin_size - bin_size 
     * @returns - array - read depth statistics array
     */
    async get_rd_stat(bin_size){
    
        let rd_stat_signal =  `rd_stat_${bin_size}_auto`
        let rd_stat;
        if (this.pytorKeys.includes(rd_stat_signal)){
            const rd_stat_dataset = await this.h5_obj.get(rd_stat_signal)
            rd_stat = await rd_stat_dataset.value
        }
        return rd_stat
    }

    
    async get_chr_signal(chrom, bin_size, signal_name, rd_stat){
        /* return a list of dictionary for a chromosome */
        let chr_wig = [];
        
        if (this.pytorKeys.includes(signal_name)){
            const chrom_dataset = await this.h5_obj.get(signal_name)
            
            let chrom_data = await chrom_dataset.value
            chrom_data.forEach((bin_value, bin_idx) => {
                chr_wig.push({chr:chrom, start: bin_idx*bin_size, end: (bin_idx+1) * bin_size, value: (bin_value/rd_stat[4]) *2})
            });
        }
        return chr_wig
    }


    async getBafSignals(chrom, binSize, signalName, scalingFactor = -1) {
        const chrWig1 = [];
        const chrWig2 = [];
        
        if (this.pytorKeys.includes(signalName)) {
            const chromDataset = await this.h5_obj.get(signalName);
            const chromData = await chromDataset.to_array();

            chromData.forEach((lh, binIdx) => {
                if (!isNaN(lh)) {
                    chrWig1.push({
                        chr: chrom,
                        start: binIdx * binSize,
                        end: (binIdx + 1) * binSize,
                        value: scalingFactor * (0.5 - lh)
                    });
                    if (lh !== 0.5) {
                        chrWig2.push({
                            chr: chrom,
                            start: binIdx * binSize,
                            end: (binIdx + 1) * binSize,
                            value: scalingFactor * (0.5 + lh)
                        });
                    }
                }
            });
        }
        return [chrWig1, chrWig2];
    }

}

class ParseSignals {
    /**
     * @param {string[]} signals - List of keys in pytor files.
     */
    constructor(signals) {
        this.signals = signals;
    }

    getAllBins() {
        const rdBins = this.getRdBins();
        const snpBins = this.getSnpBins();
        return [...new Set([...rdBins, ...snpBins])].sort((a, b) => a - b);;
    }

    getRdBins() {
        return this.extractBins(/^his_rd_p_(.*)_(\d+)$/);
    }

    getSnpBins() {
        return this.extractBins(/^snp_likelihood_(.*)_(\d+)_mask$/);
    }

    extractBins(regex) {
        return [...new Set(
            this.signals
                .map(val => val.match(regex))
                .filter(match => match !== null)
                .map(match => Number(match[2]))
        )];
    }
}

function fixString(strings) {

    return strings.map(s => s.substr(0,s.indexOf('\0')))

}

// function to_array(value, shape) {
//     const { json_value, metadata } = this;
//     const { shape } = metadata;
//     if (!isIterable(json_value) || typeof json_value === "string") {
//         return json_value;
//     }
//     let nested = create_nested_array(json_value, shape);
//     return nested;
// }

function create_nested_array(value, shape) {
    // check that shapes match:
    const total_length = value.length;
    const dims_product = shape.reduce((previous, current) => (previous * current), 1);
    if (total_length !== dims_product) {
        console.warn(`shape product: ${dims_product} does not match length of flattened array: ${total_length}`);
    }
    // Get reshaped output:
    let output = value;
    const subdims = shape.slice(1).reverse();
    for (let dim of subdims) {
        // in each pass, replace input with array of slices of input
        const new_output = [];
        const { length } = output;
        let cursor = 0;
        while (cursor < length) {
            new_output.push(output.slice(cursor, cursor += dim));
        }
        output = new_output;
    }
    return output;
}



export default HDF5Reader
