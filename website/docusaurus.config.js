// @ts-check
import {themes as prismThemes} from 'prism-react-renderer';
import remarkMath from 'remark-math';
import rehypeKatex from 'rehype-katex';

/** @type {import('@docusaurus/types').Config} */
const config = {
  title: 'pciSeq',
  tagline: 'Probabilistic cell typing by in situ sequencing',
  favicon: 'img/favicon.ico',
    
    // --- Where the site is served ---------------------------------------
  // Project site at https://acycliq.github.io/pciSeq_3d/
  url: 'https://acycliq.github.io',
  baseUrl: '/pciSeq_3d/',

  organizationName: 'acycliq',
  projectName: 'pciSeq_3d',

  onBrokenLinks: 'warn',
  markdown: {
    mermaid: true,
    hooks: {onBrokenMarkdownLinks: 'warn'},
  },

  i18n: {defaultLocale: 'en', locales: ['en']},

  themes: ['@docusaurus/theme-mermaid'],

  presets: [
    [
      'classic',
      /** @type {import('@docusaurus/preset-classic').Options} */
      ({
        docs: {
          routeBasePath: '/',          // docs ARE the site (served at root)
          sidebarPath: './sidebars.js',
          editUrl: 'https://github.com/acycliq/pciSeq_3d/tree/dev_3d/website/',
          remarkPlugins: [remarkMath],
          rehypePlugins: [rehypeKatex],
        },
        blog: false,
        theme: {customCss: './src/css/custom.css'},
      }),
    ],
  ],

  // KaTeX stylesheet for the maths.
  stylesheets: [
    {
      href: 'https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css',
      type: 'text/css',
      integrity:
        'sha384-nB0miv6/jRmo5UMMR1wu3Gz6NLsoTkbqJghGIsx//Rlm+ZU03BU6SQNC66uf4l5+',
      crossorigin: 'anonymous',
    },
  ],

  themeConfig:
    /** @type {import('@docusaurus/preset-classic').ThemeConfig} */
    ({
      colorMode: {defaultMode: 'dark', respectPrefersColorScheme: false},
      announcementBar: {
        id: 'under_review',
        content:
          '🚧 These docs are a work in progress and the text is still under review. Expect changes.',
        backgroundColor: '#b45309',
        textColor: '#ffffff',
        isCloseable: false,
      },
      navbar: {
        title: 'pciSeq',
        items: [
          {type: 'docSidebar', sidebarId: 'docsSidebar', position: 'left', label: 'Docs'},
          {
            href: 'https://github.com/acycliq/pciSeq_3d',
            label: 'GitHub',
            position: 'right',
          },
        ],
      },
      footer: {
        style: 'dark',
        links: [
          {
            title: 'Docs',
            items: [
              {label: 'Introduction', to: '/'},
              {label: 'How it works', to: '/how-it-works/overview'},
            ],
          },
          {
            title: 'Project',
            items: [
              {label: 'pciSeq', href: 'https://github.com/acycliq/pciSeq_3d'},
              {label: 'pciSeq Viewer', href: 'https://github.com/acycliq/pciSeq_viewer'},
            ],
          },
        ],
        copyright: `Copyright © ${new Date().getFullYear()} acycliq. Built with Docusaurus.`,
      },
      prism: {
        theme: prismThemes.github,
        darkTheme: prismThemes.dracula,
        additionalLanguages: ['python', 'bash', 'json'],
      },
    }),
};

export default config;