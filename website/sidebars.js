// @ts-check

/**
 * Sidebar for the pciSeq documentation.
 * Ordered to follow how someone would learn the package:
 * what it is -> how the algorithm works, block by block -> the formal model -> reference.
 *
 * @type {import('@docusaurus/plugin-content-docs').SidebarsConfig}
 */
const sidebars = {
  docsSidebar: [
    'intro',
    {
      type: 'category',
      label: 'How it works',
      link: {type: 'doc', id: 'how-it-works/overview'},
      items: [
        'how-it-works/misread-density',
        'how-it-works/warping-the-reference',
        'how-it-works/cell-to-celltype',
        'how-it-works/spots-to-cells',
      ],
    },
    {
      type: 'category',
      label: 'The model',
      link: {type: 'doc', id: 'the-model/overview'},
      items: [
        'the-model/misread-density',
        'the-model/cell-scale-theta',
        'the-model/scale-factors',
        'the-model/cell-class',
        'the-model/spot-assignment',
        'the-model/gene-gene',
        'the-model/errata',
        'the-model/appendix-self-consistency',
      ],
    },
  ],
};

export default sidebars;