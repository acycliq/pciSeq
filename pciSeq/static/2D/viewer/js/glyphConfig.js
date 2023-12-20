// Set here the color scheme and the marker(glyph) shape for the gene panel you are working
// with. These gene-specific shapes come into action at the deep zoom levels, when Leaflet.js
// takes over.
//
// There is switch for that called zoomSwitch, see line 39 of index.js together with the moveend
// callback, see lines 459 and 511 of dapi.js.
//
// Outside this deep zoom zone all genes have the same shape, a solid circle, see
// line 11, stage_markers_patched.js
//
// Note that leaflet doesnt get on very well with lots and lots of datapoints, if you notice
// some slowness on the viewer maybe changing the zoomSwitch from its default value of 7 to 8
// will probably help.
//
// Only the following shapes are supported:
//   star6
//   star5
//   diamond
//   square
//   triangleUp
//   triangleDown
//   triangleRight
//   triangleLeft
//   cross
//   plus
//   asterisk
//   circle
//   point
//
// These are designed in glyphPaths.js and introduced to the viewer from glyphs.js (no need to edit them
// unless you want to do an extra shape).


function glyphSettings()
{
    var out = [

         {gene: 'Snca',         color: '#0000FF',   glyphName: 'plus'},
         {gene: 'Cplx2',        color: '#0000FF',   glyphName: 'point'},
         {gene: 'Lhx6',         color: '#0000FF',   glyphName: 'square'},
         {gene: 'Col25a1',      color: '#0000FF',   glyphName: 'triangleUp'},
         {gene: 'Pnoc',         color: '#0000FF',   glyphName: 'triangleRight'},
         {gene: 'Rab3c',        color: '#0000FF',   glyphName: 'triangleLeft'},
         {gene: 'Gad1',         color: '#0000FF',   glyphName: 'star5'},
         {gene: 'Slc6a1',       color: '#0000FF',   glyphName: 'star6'},
         {gene: 'Th',           color: '#00B3FF',   glyphName: 'plus'},
         {gene: 'Crhbp',        color: '#00B3FF',   glyphName: 'circle'},
         {gene: 'Sst',          color: '#00B3FF',   glyphName: 'asterisk'},
         {gene: 'Npy',          color: '#00B3FF',   glyphName: 'point'},
         {gene: 'Synpr',        color: '#00B3FF',   glyphName: 'cross'},
         {gene: 'Chodl',        color: '#00B3FF',   glyphName: 'square'},
         {gene: 'Cort',         color: '#00B3FF',   glyphName: 'diamond'},
         {gene: 'Reln',         color: '#00B3FF',   glyphName: 'triangleUp'},
         {gene: 'Serpini1',     color: '#00B3FF',   glyphName: 'triangleLeft'},
         {gene: 'Satb1',        color: '#00B3FF',   glyphName: 'triangleRight'},
         {gene: 'Grin3a',       color: '#00B3FF',   glyphName: 'star5'},
         {gene: 'Tac1',         color: '#5C33FF',   glyphName: 'circle'},
         {gene: 'Pvalb',        color: '#5C33FF',   glyphName: 'asterisk'},
         {gene: 'Kcnip2',       color: '#5C33FF',   glyphName: 'square'},
         {gene: 'Thsd7a',       color: '#5C33FF',   glyphName: 'diamond'},
         {gene: 'Cox6a2',       color: '#5C33FF',   glyphName: 'triangleDown'},
         {gene: 'Chrm2',        color: '#5C33FF',   glyphName: 'star5'},
         {gene: 'Id2',          color: '#FF00E6',   glyphName: 'plus'},
         {gene: 'Hapln1',       color: '#FF00E6',   glyphName: 'circle'},
         {gene: 'Gabrd',        color: '#FF00E6',   glyphName: 'asterisk'},
         {gene: 'Cryab',        color: '#FF00E6',   glyphName: 'cross'},
         {gene: 'Kit',          color: '#FF00E6',   glyphName: 'square'},
         {gene: 'Ndnf',         color: '#FF00E6',   glyphName: 'diamond'},
         {gene: 'Nos1',         color: '#FF00E6',   glyphName: 'triangleUp'},
         {gene: 'Lamp5',        color: '#FF00E6',   glyphName: 'triangleRight'},
         {gene: 'Cplx3',        color: '#FF00E6',   glyphName: 'star6'},
         {gene: 'Cadps2',       color: '#995C00',   glyphName: 'circle'},
         {gene: 'Cxcl14',       color: '#995C00',   glyphName: 'asterisk'},
         {gene: 'Ntng1',        color: '#995C00',   glyphName: 'square'},
         {gene: 'Cpne5',        color: '#995C00',   glyphName: 'diamond'},
         {gene: 'Rgs12',        color: '#995C00',   glyphName: 'star6'},
         {gene: 'Sncg',         color: '#FF0000',   glyphName: 'circle'},
         {gene: 'Cnr1',         color: '#FF0000',   glyphName: 'asterisk'},
         {gene: 'Cck',          color: '#FF0000',   glyphName: 'point'},
         {gene: 'Trp53i11',     color: '#FF0000',   glyphName: 'cross'},
         {gene: 'Sema3c',       color: '#FF0000',   glyphName: 'square'},
         {gene: 'Syt6',         color: '#FF0000',   glyphName: 'triangleUp'},
         {gene: 'Yjefn3',       color: '#FF0000',   glyphName: 'triangleDown'},
         {gene: 'Rgs10',        color: '#FF0000',   glyphName: 'triangleRight'},
         {gene: 'Nov',          color: '#FF0000',   glyphName: 'triangleLeft'},
         {gene: 'Kctd12',       color: '#FF0000',   glyphName: 'star5'},
         {gene: 'Slc17a8',      color: '#FF0000',   glyphName: 'star6'},
         {gene: 'Tac2',         color: '#FFC700',   glyphName: 'plus'},
         {gene: 'Npy2r',        color: '#FFC700',   glyphName: 'circle'},
         {gene: 'Calb2',        color: '#FFC700',   glyphName: 'asterisk'},
         {gene: 'Htr3a',        color: '#FFC700',   glyphName: 'point'},
         {gene: 'Slc5a7',       color: '#FFC700',   glyphName: 'cross'},
         {gene: 'Penk',         color: '#FFC700',   glyphName: 'square'},
         {gene: 'Pthlh',        color: '#FFC700',   glyphName: 'triangleUp'},
         {gene: 'Vip',          color: '#FFC700',   glyphName: 'triangleDown'},
         {gene: 'Crh',          color: '#FFC700',   glyphName: 'triangleRight'},
         {gene: 'Qrfpr',        color: '#FFC700',   glyphName: 'star5'},
         {gene: 'Zcchc12',      color: '#96B38F',   glyphName: 'plus'},
         {gene: 'Calb1',        color: '#96B38F',   glyphName: 'asterisk'},
         {gene: 'Vsnl1',        color: '#96B38F',   glyphName: 'point'},
         {gene: 'Tmsb10',       color: '#96B38F',   glyphName: 'diamond'},
         {gene: 'Rbp4',         color: '#96B38F',   glyphName: 'triangleDown'},
         {gene: 'Fxyd6',        color: '#96B38F',   glyphName: 'triangleUp'},
         {gene: '6330403K07Rik', color: '#96B38F',  glyphName: 'triangleLeft'},
         {gene: 'Scg2',         color: '#96B38F',   glyphName: 'triangleRight'},
         {gene: 'Gap43',        color: '#96B38F',   glyphName: 'star5'},
         {gene: 'Nrsn1',        color: '#96B38F',   glyphName: 'star6'},
         {gene: 'Gda',          color: '#407F59',   glyphName: 'plus'},
         {gene: 'Bcl11b',       color: '#407F59',   glyphName: 'circle'},
         {gene: 'Rgs4',         color: '#407F59',   glyphName: 'asterisk'},
         {gene: 'Slc24a2',      color: '#407F59',   glyphName: 'point'},
         {gene: 'Lphn2',        color: '#407F59',   glyphName: 'cross'},
         {gene: 'Adgrl2',       color: '#407F59',   glyphName: 'cross'},
         {gene: 'Map2',         color: '#407F59',   glyphName: 'square'},
         {gene: 'Prkca',        color: '#407F59',   glyphName: 'diamond'},
         {gene: 'Cdh13',        color: '#407F59',   glyphName: 'triangleUp'},
         {gene: 'Atp1b1',       color: '#407F59',   glyphName: 'triangleDown'},
         {gene: 'Pde1a',        color: '#407F59',   glyphName: 'triangleLeft'},
         {gene: 'Calm2',        color: '#407F59',   glyphName: 'triangleRight'},
         {gene: 'Sema3e',       color: '#407F59',   glyphName: 'star6'},
         {gene: 'Nrn1',         color: '#00FF00',   glyphName: 'asterisk'},
         {gene: 'Pcp4',         color: '#00FF00',   glyphName: 'point'},
         {gene: 'Rprm',         color: '#00FF00',   glyphName: 'plus'},
         {gene: 'Enpp2',        color: '#00FF00',   glyphName: 'cross'},
         {gene: 'Rorb',         color: '#00FF00',   glyphName: 'circle'},
         {gene: 'Rasgrf2',      color: '#00FF00',   glyphName: 'square'},
         {gene: 'Wfs1',         color: '#00FF00',   glyphName: 'diamond'},
         {gene: 'Fos',          color: '#00FF00',   glyphName: 'triangleRight'},
         {gene: 'Plcxd2',       color: '#00FF00',   glyphName: 'triangleDown'},
         {gene: 'Crym',         color: '#00FF00',   glyphName: 'triangleLeft'},
         {gene: '3110035E14Rik', color: '#00FF00',  glyphName: 'triangleUp'},
         {gene: 'Foxp2',        color: '#00FF00',   glyphName: 'star5'},
         {gene: 'Pvrl3',        color: '#00FF00',   glyphName: 'star6'},
         {gene: 'Neurod6',      color: '#44B300',   glyphName: 'plus'},
         {gene: 'Nr4a2',        color: '#44B300',   glyphName: 'circle'},
         {gene: 'Cux2',         color: '#44B300',   glyphName: 'asterisk'},
         {gene: 'Kcnk2',        color: '#44B300',   glyphName: 'point'},
         {gene: 'Arpp21',       color: '#44B300',   glyphName: 'square'},
         {gene: 'Enc1',         color: '#44B300',   glyphName: 'triangleDown'},
         {gene: 'Fam19a1',      color: '#44B300',   glyphName: 'triangleRight'},
         {gene: 'Vim',          color: '#FFFFFF',   glyphName: 'asterisk'},
         {gene: 'Slc1a2',       color: '#FFFFFF',   glyphName: 'point'},
         {gene: 'Pax6',         color: '#FFFFFF',   glyphName: 'square'},
         {gene: 'Plp1',         color: '#FFFFFF',   glyphName: 'cross'},
         {gene: 'Mal',          color: '#FFFFFF',   glyphName: 'plus'},
         {gene: 'Aldoc',        color: '#FFFFFF',   glyphName: 'circle'},
         {gene: 'Actb',         color: '#FFFFFF',   glyphName: 'triangleDown'},
         {gene: 'Sulf2',        color: '#FFFFFF',   glyphName: 'star5'},


        // Do not remove this. If a gene is missing from the settings above, use these as default
        {gene: 'Generic',       color: '#0000FF',   glyphName: 'circle'},

        ];

    return out
}
