import TrackBase from "../trackBase.js"
import HDF5IndexedReader from "./HDF5IndexedReader.js"
import {CNVpytorVCF} from "./cnvpytorVCF.js"
import FeatureSource from '../feature/featureSource.js'
import {createCheckbox} from "../igv-icons.js"
import IGVGraphics from "../igv-canvas.js"
import VariantTrack from "../variant/variantTrack.js"

// testing for highlighting CNV calls
import * as DOMUtils from "../ui/utils/dom-utils.js"
import { FileUtils} from "../../node_modules/igv-utils/src/index.js"
import TrackROISet from "../roi/trackROISet.js"

class CNVPytorTrack extends TrackBase {

    static DEFAULT_TRACK_HEIGHT = 250

    constructor(config, browser) {
        super(config, browser)
    }

    init(config) {

        // NOTE -- don't use the "defaults" convention for this track, it will not work with VariantTrack.convertToPytor()
        this.featureType = 'numeric'
        this.type = "cnvpytor"
        if (!config.max) {
            this.defaultScale = true
            this.autoscale = false
        }
        if(!config.height) config.height = CNVPytorTrack.DEFAULT_TRACK_HEIGHT

        this.type = "cnvpytor"
        this.graphType = config.graphType || "points"
        this.bin_size = config.bin_size || 100000
        this.signal_name = config.signal_name || "rd_snp"
        this.cnv_caller = config.cnv_caller || '2D'
        this.colors = config.colors || ['gray', 'black', 'green', 'blue']
        this.hasChroms = {}
        this.highlightCNV = true
    
        super.init(config)

    }

    get supportsWholeGenome() {
        return true
    }

    get_signals() {
        let signals = []

        if (this.signal_name == 'rd_snp') {
            signals = ["RD_Raw", "RD_Raw_gc_corrected", this.cnv_caller, "BAF1", "BAF2"]

        } else if (this.signal_name == 'rd') {
            signals = ["RD_Raw", "RD_Raw_gc_corrected", this.cnv_caller]

        } else if (this.signal_name == 'snp') {
            signals = ["BAF1", "BAF2"]

        } else if (this.signal_name == 'cnh') {
            signals = [this.cnv_caller]
        }
        return signals
    }

    get_signal_colors() {

        let signal_colors = [
            {singal_name: 'RD_Raw', color: this.colors[0]},
            {singal_name: 'RD_Raw_gc_coor', color: this.colors[1]},
            {singal_name: 'ReadDepth', color: this.colors[2]},
            {singal_name: '2D', color: this.colors[2]},
            {singal_name: 'BAF1', color: this.colors[3]},
            {singal_name: 'BAF2', color: this.colors[3]},
        ]
        return signal_colors
    }

    async postInit() {

        if (this.config.format == 'vcf') {

            let allVariants
            if (this.featureSource) {
                allVariants = Object.values(this.featureSource.getAllFeatures()).flat()
            } else {
                this.featureSource = this.featureSource || FeatureSource(this.config, this.browser.genome)
                this.header = await this.getHeader()
                allVariants = this.featureSource.reader.features
            }

            const refGenome = this.browser.config.genome
            
            // Initializing CNVpytorVCF class
            const cnvpytor_obj = new CNVpytorVCF(allVariants, this.bin_size, refGenome)
            let wigFeatures
            let bafFeatures
            this.wigFeatures_obj = {}
            this.wigFeatures_obj[this.bin_size] = {}

            let dataWigs

            if (this.cnv_caller == '2D') {

                dataWigs = await cnvpytor_obj.read_rd_baf('2D')

                wigFeatures = dataWigs[0]
                bafFeatures = dataWigs[1]
                this.wigFeatures_obj[this.bin_size]['2D'] = wigFeatures[2]

                this.available_callers = ['2D']
            } else {
                dataWigs = await cnvpytor_obj.read_rd_baf()
                wigFeatures = dataWigs[0]
                bafFeatures = dataWigs[1]
                this.wigFeatures_obj[this.bin_size]['ReadDepth'] = wigFeatures[2]
                this.available_callers = ['ReadDepth']
            }

            this.wigFeatures_obj[this.bin_size]['RD_Raw'] = wigFeatures[0]
            this.wigFeatures_obj[this.bin_size]['RD_Raw_gc_coor'] = wigFeatures[1]
            this.wigFeatures_obj[this.bin_size]['BAF1'] = bafFeatures[0]
            this.wigFeatures_obj[this.bin_size]['BAF2'] = bafFeatures[1]

            this.available_bins = [this.bin_size]

            this.set_available_callers()

        } else {
            this.cnvpytor_obj = new HDF5IndexedReader(this.config.url, this.bin_size)
            // get chrom list that currently user viewing
            let chroms = [ ...new Set(this.browser.referenceFrameList.map(val => val.chr))]
            
            let aliasRecord = this.getAliasChromsList(chroms)
            // console.log("aliasRecord", aliasRecord)
            // console.log("chrom", chroms)
            //this.wigFeatures_obj = await this.cnvpytor_obj.get_rd_signal(this.bin_size, this.cnv_caller, aliasRecord)
            
            await this.cnvpytor_obj.get_rd_signal(this.bin_size, this.cnv_caller, aliasRecord)
            this.wigFeatures_obj = this.cnvpytor_obj.allWigFeatures
            
            // Save the processed chroms names to check later for the availability
            this.update_hasChroms(this.wigFeatures_obj, chroms)

            // console.log("this.hasChroms : ", this.hasChroms)
            this.available_bins = this.cnvpytor_obj.availableBins
            // reset the bin size if its not exits
            if(! this.available_bins.includes(this.bin_size)){
                this.bin_size = this.available_bins.at(-1);    
            }

            this.available_callers = this.cnvpytor_obj.callers
            this.set_available_callers()
        }

        this.tracks = []
        const p = []

        this.signals = this.get_signals()
        this.signal_colors = this.get_signal_colors()

        for (let bin_size in this.wigFeatures_obj) {
            let i = 0
            for (const [signal_name, wig] of Object.entries(this.wigFeatures_obj[bin_size])) {
                
                if (!Array.isArray(wig)) {continue}
        
                if (this.signals.includes(signal_name)) {
                    let tconf = {}
                    tconf.type = "wig"
                    tconf.isMergedTrack = true
                    tconf.features = wig
                    tconf.name = signal_name
                    tconf.color = this.signal_colors.filter(x => x.singal_name === signal_name).map(x => x.color)
                    const t = await this.browser.createTrack(tconf)
                    if (t) {
                        t.autoscale = false     // Scaling done from merged track
                        this.tracks.push(t)
                    } else {
                        console.warn("Could not create track " + tconf)
                    }

                    if (typeof t.postInit === 'function') {
                        p.push(t.postInit())
                    }
                    i++
                }
            }

        }

        this.flipAxis = this.config.flipAxis ? this.config.flipAxis : false
        this.logScale = this.config.logScale ? this.config.logScale : false
        this.autoscale = this.config.autoscale
        if (!this.autoscale) {
            this.dataRange = {
                min: this.config.min || 0,
                max: this.config.max
            }
        }
        for (let t of this.tracks) {
            t.autoscale = false
            t.dataRange = this.dataRange
        }
        
        this.default_filter_criteria = {p_range: [0, 1], q0_range: [0, 1], dG:0, cf:0.00, baf: [0.01, 1], bins: 4} 
        // for meanshift 
        this.filter_criteria = {p_range: [0, 1], q0_range: [0, 1] , dG:0}
        // 2d callers
        this.filter_criteria = {p_range: [0, 1], q0_range: [0, 1], cf:0.00, baf: [0.01, 1], bins: 1}
        this.get_roi_features()
        
        return Promise.all(p)
    }

    get_roi_features() {
        // Initialize roiSets
        this.roiSets = [];
    
        // Fetch CNV calls
        const cnvCalls = this.cnvpytor_obj.CNVcalls?.[this.bin_size]?.[this.cnv_caller];
    
        if (this.highlightCNV && cnvCalls) {
            const filteredCNV_Calls = this.get_filtered_calls(cnvCalls);
            if (filteredCNV_Calls.length > 0) {
                this.roiSets = filteredCNV_Calls.map(r => new TrackROISet(r, this.browser.genome));
            }
        }
    
        // Add additional ROIs from config if available
        if (this.config.roi?.length) {
            this.roiSets.push(...this.config.roi.map(r => new TrackROISet(r, this.browser.genome)));
        }
    }
    
    get_filtered_calls(cnvCalls){
        let filteredCNV_Calls = [];
        const colorMap = {
            deletion: "rgba(250,0,0,0.25)",
            duplication: "rgba(3,52,249,0.25)",
            cnnloh: "rgba(3,250,0,0.25)"
        };
        const isInRange = (value, min, max) => min <= value && value <= max;

        for (const [cnv_type, cnv_rows] of Object.entries(cnvCalls)) {
            
            // apply the filtering here
            let filteredCNV_Rows = cnv_rows.filter(row => 
                isInRange(row.p_val, this.filter_criteria.p_range[0], this.filter_criteria.p_range[1]) &&
                isInRange(row.Q0, this.filter_criteria.q0_range[0], this.filter_criteria.q0_range[1]) &&
        
                Math.abs(row.cnv - 1) * 2 > this.filter_criteria.cf &&
                (row.baf === undefined || isInRange(row.baf, this.filter_criteria.baf[0], this.filter_criteria.baf[1]) ) &&
                (row.bins === undefined || row.bins >= this.filter_criteria.bins)
            )
            
            if (filteredCNV_Rows.length > 0) {
                
                filteredCNV_Calls.push({
                    name: cnv_type,
                    features: filteredCNV_Rows,
                    color: colorMap[cnv_type] || "rgba(0,0,0,0.25)" // Default color (optional)
                });
            }
        }
        
        return filteredCNV_Calls
    }

    
    getAliasChromsList(chroms){
        let aliasRecord = chroms.map(chr => {
            let records = this.browser.genome.chromAlias.aliasRecordCache.get(chr)
            return Object.values(records)
        })
        aliasRecord = aliasRecord.flat()
        return aliasRecord
    }
    
    set_available_callers() {
        if (!this.available_callers.includes(this.cnv_caller)) {
            if (this.available_callers.length > 0) {
                this.cnv_caller = this.available_callers[0]
            } else {
                this.cnv_caller = null
            }
        }
    }

    async getHeader() {

        if (!this.header) {
            if (typeof this.featureSource.getHeader === "function") {
                const header = await this.featureSource.getHeader()
                if (header) {
                    this.callSets = header.callSets || []
                }
                this.header = header
            }
            this.sampleKeys = this.callSets ? this.callSets.map(cs => cs.sample) : []
            this.sampleNames = this.sampleKeys
        }

        return this.header
    }

    get height() {
        return this._height
    }

    set height(h) {
        this._height = h
        if (this.tracks) {
            for (let t of this.tracks) {
                t.height = h
                t.config.height = h
            }
        }
    }

    menuItemList() {
        let items = []

        if (this.flipAxis !== undefined) {
            items.push({
                label: "Flip y-axis",
                click: function flipYAxisHandler() {
                    this.flipAxis = !this.flipAxis
                    this.trackView.repaintViews()
                }
            })
        }

        items = items.concat(this.numericDataMenuItems())
        items.push('<hr/>')
        items.push("Highlight CNV Calls")
        for (let hc of [true, false]){
            items.push({
                element: createCheckbox(hc, hc === this.highlightCNV),
                click: async function roiHandler() {
                    if (hc !== this.highlightCNV) {
                        
                        this.highlightCNV = hc
                        this.get_roi_features()
                       
                        this.clearCachedFeatures()
                        
                        this.trackView.updateViews()
                        this.trackView.repaintViews()
                    }
                    
                }
            })
        }
        // items.push('<hr/>')
        
        items.push(this.PvalueAdjustmentMenuItem())
        items.push(this.q0AdjustmentMenuItem())
        items.push(this.BinsAdjustmentMenuItem())
        items.push(this.GapAdjustmentMenuItem())
        items.push(this.CellFreqAdjustmentMenuItem() )
        items.push(this.BafAdjustmentMenuItem())

        items.push('<hr/>')
        items.push("Bin Sizes")
        for (let rd_bin of this.available_bins) {

            items.push({
                element: createCheckbox(rd_bin, rd_bin === this.bin_size),
                click: async function binSizesHandler() {
                    if (rd_bin !== this.bin_size) {
                        this.bin_size = rd_bin
                        // data loader image
                        this.trackView.startSpinner()

                        await this.recreate_tracks(rd_bin)
                        this.clearCachedFeatures()
                        this.trackView.updateViews()
                        this.trackView.repaintViews()
                    }
                }
            })
        }
        items.push('<hr/>')
        items.push("Signal Type")

        let signal_dct = {"rd_snp": "RD and BAF Likelihood", "rd": "RD Signal", "snp": "BAF Likelihood"}
        for (let signal_name in signal_dct) {

            items.push({
                element: createCheckbox(signal_dct[signal_name], signal_name === this.signal_name),
                click: async function signalTypeHandler() {
                    this.signal_name = signal_name
                    await this.recreate_tracks(this.bin_size)
                    this.clearCachedFeatures()
                    this.trackView.updateViews()
                    this.trackView.repaintViews()

                }
            })
        }

        // cnv caller setting
        items.push('<hr/>')
        items.push("CNV caller")
        for (let cnv_caller of this.available_callers) {

            items.push({
                element: createCheckbox(cnv_caller, cnv_caller === this.cnv_caller),
                click: async function cnvCallerHandler() {
                    if (cnv_caller !== this.cnv_caller) {
                        this.cnv_caller = cnv_caller
                        // data loader image
                        this.trackView.startSpinner()

                        await this.recreate_tracks(this.bin_size)
                        this.clearCachedFeatures()
                        this.trackView.updateViews()
                        this.trackView.repaintViews()
                    }
                }
            })
        }

        // download calls
        items.push('<hr/>')
        items.push("Download CNV calls")
        items.push({
            label: "\u00A0\u00A0 All calls",
            click: function downloadCallHandler() {

                const path = 'cnvpytortrack_igvjs.tsv';
                let CNV_calls = this.cnvpytor_obj.CNVcalls[this.bin_size][this.cnv_caller];
                this.getTSVformattedCalls(path, CNV_calls)
            }
        })
        items.push({
            label: "\u00A0\u00A0 Filtered calls",
            click: function downloadCallHandler() {

                const path = 'cnvpytortrack_filtered_igvjs.tsv';
                let CNV_calls = this.cnvpytor_obj.CNVcalls[this.bin_size][this.cnv_caller];


                let filtered_roi = this.get_filtered_calls(CNV_calls)                
                if (filtered_roi){
                    let download_calls = {}
                    for (let call of filtered_roi){
                        if(call.features.length > 0){
                            download_calls[call.name] = call.features
                        }
                    }
                    
                    this.getTSVformattedCalls(path, download_calls)
                }
            }
        })

        // variant track conversion -- only if track was originally created from a VariantTrack
        if (this.variantState) {
            items.push('<hr/>')
            for (let cnv_caller of this.available_callers) {
                items.push({
                    label: 'Convert to variant track',
                    click: () => {
                        this.convertToVariant()
                    }
                })

            }
        }


        return items
    }
    getTSVformattedCalls(path, calls){
        // Flattening the JSON structure into an array of rows
        const rows = [];
        const headers = new Set(); // Collect unique headers dynamically

        // Iterate over each category (deletion, duplication)
        for (const [category, entries] of Object.entries(calls)) {
            entries.forEach(entry => {
                entry["type"] = category; // Add a column to differentiate types
                rows.push(entry);
                Object.keys(entry).forEach(key => headers.add(key)); // Collect all headers
            });
        }

        // Define the required column order
        const requiredOrder = ["type", "chr", "start", "end", "cnv", "size"];

        // Collect remaining headers dynamically
        const otherHeaders = [...headers].filter(h => !requiredOrder.includes(h));

        // Final column order: required columns first, then remaining columns
        const finalHeaders = [...requiredOrder, ...otherHeaders];

        // Convert data to TSV format
        const tsvData = [
            finalHeaders.join("\t"), // Header row
            ...rows.map(row => finalHeaders.map(field => row[field] ?? "").join("\t")) // Data rows
        ].join("\n");

        // Create a Blob and Object URL
        const blob = new Blob([tsvData], { type: "text/tab-separated-values" });
        const dataUrl = URL.createObjectURL(blob);

        // Trigger the download
        FileUtils.download(path, dataUrl);

        // Revoke the URL after a short delay to release memory
        setTimeout(() => URL.revokeObjectURL(dataUrl), 1000);
    }

    PvalueAdjustmentMenuItem() {
        const element = DOMUtils.div()
        element.innerText = 'Set Max P value'

        function dialogPresentationHandler(e) {
            const callback = alpha => {
            
                this.filter_criteria.p_range[1] = Math.max(0.0000, alpha)
                this.get_roi_features()

                // this.clearCachedFeatures()
                this.trackView.updateViews()
                this.trackView.repaintViews()
            }

            const config ={
                label: 'Max P value',
                value: this.filter_criteria.p_range[1],
                min: 0.00,
                max: this.default_filter_criteria.p_range[1],
                precision: 4,
                scaleFactor: 1000,
                callback
            }

            this.browser.sliderDialog.present(config, e)
            
        }
        
        return {element, dialog: dialogPresentationHandler}
    }

    q0AdjustmentMenuItem() {
        const element = DOMUtils.div()
        element.innerText = 'Set Max Q0 value'

        function dialogPresentationHandler(e) {
            const callback = alpha => {
            
                this.filter_criteria.q0_range[1] = Math.max(0.0000, alpha)

                this.get_roi_features()

                // this.clearCachedFeatures()
                this.trackView.updateViews()
                this.trackView.repaintViews()
            }

            const config ={
                label: 'Max Q0 value',
                value: this.filter_criteria.q0_range[1],
                min: 0.00,
                max: this.default_filter_criteria.q0_range[1],
                precision: 4,
                scaleFactor: 1000,
                callback
            }

            this.browser.sliderDialog.present(config, e)
            
        }
        
        return {element, dialog: dialogPresentationHandler}
    }

    GapAdjustmentMenuItem() {
        const element = DOMUtils.div()
        element.innerText = 'Set Min Distance from GAP'

        function dialogPresentationHandler(e) {
            const callback = alpha => {
                
                this.filter_criteria.dG = Math.max(0.0000, alpha)

                this.get_roi_features()
                //this.clearCachedFeatures()
                this.trackView.updateViews()
                // this.repaintViews()
                this.trackView.repaintViews()
            }

            const config ={
                label: 'Min Distance from GAP',
                value: this.default_filter_criteria.dG,
                min: 0.00,
                max: 200000,
                precision: 0,
                scaleFactor: 1,
                callback
            }
            
            this.browser.sliderDialog.present(config, e)
            
        }
        return {element, dialog: dialogPresentationHandler}
    }

    BinsAdjustmentMenuItem() {
        const element = DOMUtils.div()
        element.innerText = 'Set Min Bins'

        function dialogPresentationHandler(e) {
            const callback = alpha => {
                
                this.filter_criteria.bins = Math.max(0.0000, alpha)

                this.get_roi_features()
                
                //this.clearCachedFeatures()
                this.trackView.updateViews()
                // this.repaintViews()
                this.trackView.repaintViews()
            }

            const config ={
                label: 'Min Bin count',
                value: this.filter_criteria.bins,
                min: 0.00,
                max: 20,
                precision: 0,
                scaleFactor: 1,
                callback
            }   
            this.browser.sliderDialog.present(config, e)
            
        }
        return {element, dialog: dialogPresentationHandler}
    }

    CellFreqAdjustmentMenuItem() {
        const element = DOMUtils.div()
        element.innerText = 'Set Min Cell frequency'

        function dialogPresentationHandler(e) {
            const callback = alpha => {
                
                this.filter_criteria.cf = Math.max(0.0000, alpha)

                this.get_roi_features()
                //this.clearCachedFeatures()
                this.trackView.updateViews()
                // this.repaintViews()
                this.trackView.repaintViews()
            }

            const config ={
                label: 'Min Cell frequency',
                value: this.filter_criteria.cf,
                min: 0.00,
                max: 1,
                precision: 0,
                scaleFactor: 100,
                callback
            }   
            this.browser.sliderDialog.present(config, e)
            
        }
        return {element, dialog: dialogPresentationHandler}
    }

    BafAdjustmentMenuItem() {
        const element = DOMUtils.div()
        element.innerText = 'Set Min BAF'

        function dialogPresentationHandler(e) {
            const callback = alpha => {
                
                this.filter_criteria.baf[0] = Math.max(0.0000, alpha)

                this.get_roi_features()
                //this.clearCachedFeatures()
                this.trackView.updateViews()
                // this.repaintViews()
                this.trackView.repaintViews()
            }

            const config = {
                label: 'Min BAF',
                value: this.filter_criteria.baf[0],
                min: 0.00,
                max: 0.5,
                precision: 4,
                scaleFactor: 1000,
                callback
            }   
            this.browser.sliderDialog.present(config, e)
            
        }
        return {element, dialog: dialogPresentationHandler}
    }
    
    async recreate_tracks(bin_size) {
        this.tracks = []
        const p = []

        // if (!(bin_size in this.wigFeatures_obj)) {
        if (!this.cnvpytor_obj.allWigFeatures[bin_size]) {
            // console.log('Update track 2', this.hasChroms)
            if(Object.keys(this.hasChroms).length > 0) {
                // console.log("test2")
                let chroms = [ ...new Set(this.browser.referenceFrameList.map(val => val.chr))]
                if(chroms[0] == "all"){
                    chroms = this.browser.genome.chromosomeNames
                }
                
                // this.wigFeatures_obj = {...this.wigFeatures_obj, ...await this.cnvpytor_obj.get_rd_signal(bin_size, this.cnv_caller, chroms)}
                
                await this.cnvpytor_obj.get_rd_signal(bin_size, this.cnv_caller, chroms)
                
                this.wigFeatures_obj = this.cnvpytor_obj.allWigFeatures
                this.update_hasChroms(this.wigFeatures_obj, chroms)
                
            } else{
                // console.log('update bin size', bin_size)
                // if (this.cnvpytor_obj.allWigFeatures[bin_size][this.callers] == undefined) {
                await this.cnvpytor_obj.get_rd_signal(bin_size, this.cnv_caller)
                // this.wigFeatures_obj = {...this.wigFeatures_obj, ...await this.cnvpytor_obj.get_rd_signal(bin_size, this.cnv_caller)}
                
            }
            
        }
        if (!this.cnvpytor_obj.allWigFeatures[bin_size][this.cnv_caller]){

            let chroms = [ ...new Set(this.browser.referenceFrameList.map(val => val.chr))]
            if(chroms[0] == "all"){
                chroms = this.browser.genome.chromosomeNames
            }
            // cnv_caller, rdChromosomes, bin_size
            await this.cnvpytor_obj.getCallerWigfeatures(this.cnv_caller, chroms, bin_size)
        }

        this.get_roi_features()

        this.signals = this.get_signals()
        this.signal_colors = this.get_signal_colors()

        let i = 0

        for (const [signal_name, wig] of Object.entries(this.wigFeatures_obj[bin_size])) {
            if (this.signals.includes(signal_name)) {
                let tconf = {}
                tconf.type = "wig"
                tconf.isMergedTrack = true
                tconf.features = wig
                tconf.name = signal_name
                tconf.color = this.signal_colors.filter(x => x.singal_name === signal_name).map(x => x.color)
                const t = await this.browser.createTrack(tconf)
                if (t) {
                    t.autoscale = false     // Scaling done from merged track
                    this.tracks.push(t)
                } else {
                    console.warn("Could not create track " + tconf)
                }

                if (typeof t.postInit === 'function') {
                    p.push(t.postInit())
                }
                i++
            }

        }

        this.flipAxis = this.config.flipAxis ? this.config.flipAxis : false
        this.logScale = this.config.logScale ? this.config.logScale : false
        this.autoscale = this.config.autoscale
        if (!this.autoscale) {
            this.dataRange = {
                min: this.config.min || 0,
                max: this.config.max
            }
        }
        for (let t of this.tracks) {
            t.autoscale = false
            t.dataRange = this.dataRange
        }
        return Promise.all(p)
    }

    update_hasChroms(wigFeatures, chroms){
        for (let binSize in wigFeatures){
            if (!this.hasChroms[binSize]) {
                this.hasChroms[binSize] = new Set();
            }
            chroms.forEach(item => this.hasChroms[binSize].add(item))

        }
        return this.hasChroms

    }

    async getFeatures(chr, bpStart, bpEnd, bpPerPixel) {
        // console.log("get features")
        if(Object.keys(this.hasChroms).length > 0) {
            
            // Need to find the current binSize
            if (this.hasChroms[this.bin_size].size != 0){
                let chroms = [ ...new Set(this.browser.referenceFrameList.map(val => val.chr))]
                if(chroms[0] == "all"){
                    chroms = this.browser.genome.chromosomeNames
                }
                let newChroms = chroms.filter(val => !this.hasChroms[this.bin_size].has(val))

                if(newChroms.length != 0){

                    let aliasRecords = this.getAliasChromsList(newChroms)
                    // update the hasChroms list
                    let tmp_wig = await this.cnvpytor_obj.get_rd_signal(this.bin_size, this.cnv_caller, aliasRecords)
                    this.update_hasChroms(tmp_wig, newChroms)

                    // here we need to update the wig
                    // this part is probaby not required; code improve required
                    
                    this.get_roi_features()

                    for (let bin_size in tmp_wig){
                        for (const [signal_name, wig] of Object.entries(tmp_wig[bin_size])) {
                            await this.wigFeatures_obj[bin_size][signal_name].push(...wig)
                        }
                    }

                    for (let wig of this.tracks){                       
                        await wig.featureSource.updateFeatures(this.wigFeatures_obj[this.bin_size][wig.name] )
                    }
                }
            }

        }

        if (this.tracks) {
            const promises = this.tracks.map((t) => t.getFeatures(chr, bpStart, bpEnd, bpPerPixel))
            return Promise.all(promises)
        } else {
            return undefined  // This can happen if a redraw is triggered before the track has initialized.
        }
    }

    // TODO: refactor to igvUtils.js
    getScaleFactor(min, max, height, logScale) {
        const scale = logScale ? height / (Math.log10(max + 1) - (min <= 0 ? 0 : Math.log10(min + 1))) : height / (max - min)
        return scale
    }

    computeYPixelValue(yValue, yScaleFactor) {
        return (this.flipAxis ? (yValue - this.dataRange.min) : (this.dataRange.max - yValue)) * yScaleFactor
    }

    computeYPixelValueInLogScale(yValue, yScaleFactor) {
        let maxValue = this.dataRange.max
        let minValue = this.dataRange.min
        if (maxValue <= 0) return 0 // TODO:
        if (minValue <= -1) minValue = 0
        minValue = (minValue <= 0) ? 0 : Math.log10(minValue + 1)
        maxValue = Math.log10(maxValue + 1)
        yValue = Math.log10(yValue + 1)
        return ((this.flipAxis ? (yValue - minValue) : (maxValue - yValue)) * yScaleFactor)
    }

    draw(options) {

        // const mergedFeatures = options.features    // Array of feature arrays, 1 for each track
        const mergedFeatures = options.features
        if (!mergedFeatures) return

        if (this.defaultScale) {
            if (this.signal_name == 'rd_snp') {
                this.dataRange = {
                    min: this.config.min || this.dataRange.min || -1,
                    max: this.config.max || this.dataRange.max || 5
                }
            } else if (this.signal_name == 'rd') {
                this.dataRange = {
                    min: this.config.min || this.dataRange.min || 0,
                    max: this.config.max || this.dataRange.max || 5
                }
            } else if (this.signal_name == 'snp') {
                this.dataRange = {
                    min: this.config.min || this.dataRange.min || -1,
                    max: this.config.max || this.dataRange.max || 0
                }
            }
        }

        if (this.autoscale) {
            this.dataRange = autoscale(options.referenceFrame.chr, mergedFeatures)
        }

        if (this.tracks) {
            for (let i = 0, len = this.tracks.length; i < len; i++) {
                const trackOptions = Object.assign({}, options)
                trackOptions.features = mergedFeatures[i]
                this.tracks[i].dataRange = this.dataRange
                this.tracks[i].flipAxis = this.flipAxis
                this.tracks[i].logScale = this.logScale
                if (this.graphType) {
                    this.tracks[i].graphType = this.graphType
                }
                this.tracks[i].draw(trackOptions)
            }
        }

        // guides lines
        const scaleFactor = this.getScaleFactor(this.dataRange.min, this.dataRange.max, options.pixelHeight, this.logScale)
        const yScale = (yValue) => this.logScale
            ? this.computeYPixelValueInLogScale(yValue, scaleFactor)
            : this.computeYPixelValue(yValue, scaleFactor)

        // Draw guidelines
        if (this.config.hasOwnProperty('guideLines')) {
            for (let line of this.config.guideLines) {
                if (line.hasOwnProperty('color') && line.hasOwnProperty('y') && line.hasOwnProperty('dotted')) {
                    let y = yScale(line.y)
                    let props = {
                        'strokeStyle': line['color'],
                        'strokeWidth': 1
                    }
                    if (line['dotted']) IGVGraphics.dashedLine(options.context, 0, y, options.pixelWidth, y, 5, props)
                    else IGVGraphics.strokeLine(options.context, 0, y, options.pixelWidth, y, props)
                }
            }
        }

        let props = {
            'strokeStyle': 'lightgray',
            'strokeWidth': 0.5
        }
        let y = yScale(2)
        IGVGraphics.dashedLine(options.context, 0, y, options.pixelWidth, y, 5, props)

    }

    paintAxis(ctx, pixelWidth, pixelHeight) {

        var x1,
            x2,
            y1,
            y2,
            a,
            b,
            reference,
            shim,
            font = {
                'font': 'normal 10px Arial',
                'textAlign': 'right',
                'strokeStyle': "black"
            }

        if (undefined === this.dataRange || undefined === this.dataRange.max || undefined === this.dataRange.min) {
            return
        }

        let flipAxis = (undefined === this.flipAxis) ? false : this.flipAxis

        IGVGraphics.fillRect(ctx, 0, 0, pixelWidth, pixelHeight, {'fillStyle': "rgb(255, 255, 255)"})

        reference = 0.95 * pixelWidth
        x1 = reference - 8
        x2 = reference

        //shim = 0.5 * 0.125;
        shim = .01
        y1 = y2 = shim * pixelHeight

        a = {x: x2, y: y1}

        // tick
        IGVGraphics.strokeLine(ctx, x1, y1, x2, y2, font)
        IGVGraphics.fillText(ctx, prettyPrint(flipAxis ? this.dataRange.min : this.dataRange.max), x1 + 4, y1 + 12, font)

        //shim = 0.25 * 0.125;
        y1 = y2 = (1.0 - shim) * pixelHeight

        b = {x: x2, y: y1}

        // tick
        IGVGraphics.strokeLine(ctx, x1, y1, x2, y2, font)
        IGVGraphics.fillText(ctx, prettyPrint(flipAxis ? this.dataRange.max : this.dataRange.min), x1 + 4, y1 - 4, font)

        IGVGraphics.strokeLine(ctx, a.x, a.y, b.x, b.y, font)

        function prettyPrint(number) {
            // if number >= 100, show whole number
            // if >= 1 show 1 significant digits
            // if <  1 show 2 significant digits

            // change the label for negative number to positive; For BAF likelihood section
            if(number < 0){
                return Math.abs(number)
            }

            if (number === 0) {
                return "0"
            } else if (Math.abs(number) >= 10) {
                return number.toFixed()
            } else if (number % 1 == 0) {
                return number.toFixed()
            } else if (Math.abs(number) >= 1) {
                return number.toFixed(1)
            } else {
                return number.toFixed(2)
            }
        }

        const scaleFactor = this.getScaleFactor(this.dataRange.min, this.dataRange.max, pixelHeight, this.logScale)
        const yScale = (yValue) => this.logScale
            ? this.computeYPixelValueInLogScale(yValue, scaleFactor)
            : this.computeYPixelValue(yValue, scaleFactor)

        const n = Math.ceil((this.dataRange.max - this.dataRange.min) / 10)
        for (let p = Math.ceil(this.dataRange.min + 1); p < Math.round(this.dataRange.max - 0.4); p += n) {
            const yp = yScale(p)
            IGVGraphics.strokeLine(ctx, 45, yp, 50, yp, font) // Offset dashes up by 2 pixel
            IGVGraphics.fillText(ctx, prettyPrint(flipAxis ? this.dataRange.max - p : p), 44, yp + 4, font)

        }

    }

    popupData(clickState, features) {

        const featuresArray = features || clickState.viewport.cachedFeatures

        if (featuresArray && featuresArray.length === this.tracks.length) {
            // Array of feature arrays, 1 for each track
            const popupData = []
            for (let i = 0; i < this.tracks.length; i++) {
                if (i > 0) popupData.push('<hr/>')
                popupData.push(`<div style=background-color:rgb(245,245,245);border-bottom-style:dashed;border-bottom-width:1px;padding-bottom:5px;padding-top:10px;font-weight:bold;font-size:larger >${this.tracks[i].name}</div>`)
                const trackPopupData = this.tracks[i].popupData(clickState, featuresArray[i])
                popupData.push(...trackPopupData)

            }
            return popupData
        }
    }

    /**
     * Applicable if this track was originally created from a VariantTrack -- attempt to convert it back
     * @returns {Promise<void>}
     */
    async convertToVariant() {

        if (this.variantState) {
            Object.setPrototypeOf(this, VariantTrack.prototype)
            this.init(this.variantState)
            await this.postInit()
            this.trackView.clearCachedFeatures()
            if(this.variantState.trackHeight) {
                this.trackView.setTrackHeight(this.variantState.trackHeight)
            }
            this.trackView.checkContentHeight()
            this.trackView.updateViews()
            delete this.variantState
        }
    }

}

function autoscale(chr, featureArrays) {

    let min = 0
    let max = -Number.MAX_VALUE
    for (let features of featureArrays) {
        for (let f of features) {
            if (typeof f.value !== 'undefined' && !Number.isNaN(f.value)) {
                min = Math.min(min, f.value)
                max = Math.max(max, f.value)
            }
        }
    }
    return {min: min, max: max}
}

export default CNVPytorTrack
