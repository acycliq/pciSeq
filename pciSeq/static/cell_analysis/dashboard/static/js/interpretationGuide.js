export class InterpretationGuide {
    constructor(svg, width, height, isReferenceGuide = false) {  // Add parameter
        this.svg = svg;
        this.width = width;
        this.height = height;
        this.padding = 5;
        this.lineHeight = 15;
        this.isReferenceGuide = isReferenceGuide;  // Add property
    }

    getText(currentUserClass, currentAssignedClass) {
        // Return different text based on plot type
        if (this.isReferenceGuide) {
            return [
                "• Points on diagonal: Perfect match between estimated and actual counts",
                "• Points above diagonal: pciSeq overestimates gene counts",
                "• Points below diagonal: pciSeq underestimates gene counts",
                // "• Distance from diagonal: Magnitude of estimation error"
            ];
        }
        // Original guide text for main plot
        return [
            `• Genes on diagonal: Contribute equally to both cell types`,
            `• Genes above diagonal: Support classification as ${currentUserClass}`,
            `• Genes below diagonal: Support classification as ${currentAssignedClass}`,
            `• Distance from diagonal: Strength of support for one type over the other`
        ];
    }

    // Keep all other methods exactly as they are
    update(currentUserClass, currentAssignedClass) {
        const guideText = this.getText(currentUserClass, currentAssignedClass);
        let guide = this.svg.select('.interpretation-guide');

        if (guide.empty()) {
            guide = this.createGuide(guideText);
        } else {
            this.updateGuideText(guide, guideText);
        }
    }

    createGuide(guideText) {
        // Keep existing implementation
        const guide = this.svg.append("g")
            .attr("class", "interpretation-guide")
            .attr("transform", `translate(${this.width - 10}, ${this.height - 10})`);

        guide.selectAll("text")
            .data(guideText)
            .enter()
            .append("text")
            .attr("x", 0)
            .attr("y", (d, i) => i * this.lineHeight)
            .style("text-anchor", "start")
            .style("font-size", "12px")
            .text(d => d);

        this.addBackgroundRect(guide);
        return guide;
    }

    addBackgroundRect(guide) {
        // Keep existing implementation
        const guideBBox = guide.node().getBBox();
        guide.insert("rect", ":first-child")
            .attr("x", guideBBox.x - this.padding)
            .attr("y", guideBBox.y - this.padding)
            .attr("width", guideBBox.width + (this.padding * 2))
            .attr("height", guideBBox.height + (this.padding * 2))
            .attr("fill", "rgba(255, 223, 186, 0.7)");

        guide.attr("transform",
            `translate(${this.width - guideBBox.width - 15}, ${this.height - guideBBox.height - 15})`);
    }

    updateGuideText(guide, guideText) {
        // Keep existing implementation
        guide.selectAll("text")
            .data(guideText)
            .text(d => d);
    }
}
