export function preparePlotData(geneNames, contrData, assignedClass, userClass) {
    // Maps gene names to their corresponding x,y coordinates
    return geneNames.map((gene, index) => ({
        name: gene,
        x: contrData[assignedClass][index],
        y: contrData[userClass][index]
    }));
}

export function prepareSelectorOptions(classes, classProbs, currentAssignedClass) {
    // Filters out the assigned class and formats the dropdown options
    return classes
        .filter(c => c !== currentAssignedClass)
        .map(c => ({
            value: c,
            label: `${c} (${(classProbs[c] * 100).toFixed(2)}%)`,
            probability: classProbs[c]
        }));
}
