// Save this as: static/js/components/InterpretationGuide.js

/**
 * Component for rendering and managing the interpretation guide overlay
 */

export class InterpretationGuide {
    constructor(svg, width, height) {
        this.svg = svg;
        this.width = width;
        this.height = height;
        this.padding = 5;
        this.lineHeight = 15;
    }

    getText(currentUserClass, currentAssignedClass) {
        // Defines the guide text with dynamic class names
        return [
            `• Genes on diagonal: Contribute equally to both cell types`,
            `• Genes above diagonal: Support classification as ${currentUserClass}`,
            `• Genes below diagonal: Support classification as ${currentAssignedClass}`,
            `• Distance from diagonal: Strength of support for one type over the other`
        ];
    }

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
        // Creates the initial guide group
        const guide = this.svg.append("g")
            .attr("class", "interpretation-guide")
            .attr("transform", `translate(${this.width - 10}, ${this.height - 10})`);

        // Add text elements
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
        // Add semi-transparent background
        const guideBBox = guide.node().getBBox();
        guide.insert("rect", ":first-child")
            .attr("x", guideBBox.x - this.padding)
            .attr("y", guideBBox.y - this.padding)
            .attr("width", guideBBox.width + (this.padding * 2))
            .attr("height", guideBBox.height + (this.padding * 2))
            .attr("fill", "rgba(255, 223, 186, 0.7)");

        // Position the guide in the bottom-right corner
        guide.attr("transform",
            `translate(${this.width - guideBBox.width - 15}, ${this.height - guideBBox.height - 15})`);
    }

    updateGuideText(guide, guideText) {
        // Updates existing guide text
        guide.selectAll("text")
            .data(guideText)
            .text(d => d);
    }
}
