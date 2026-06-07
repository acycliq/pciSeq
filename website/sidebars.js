// @ts-check

/**
 * Sidebar for the pciSeq documentation.
 * Ordered to follow how someone would learn the package:
 * what it is -> how the algorithm works, block by block -> reference.
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
  ],
};

export default sidebars;