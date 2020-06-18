// var baseTree = {
//     label: 'Base Layers',
//     children: [
//         {
//             label: 'World &#x1f5fa;',
//             children: [
//                 { label: 'OpenStreetMap', layer: osm },
//                 { label: 'Esri', layer: esri },
//                 { label: 'Google Satellite', layer: g_s },
//                 /* ... */
//             ]
//         },
//         {
//             label: 'Europe',
//             children: [
//                 { label: 'France', layer: france },
//                 { label: 'Germany', layer: germany },
//                 { label: 'Spain', layer: spain },
//                 /* ... */
//             ]
//         },
//         {
//             label: 'USA',
//             children: [
//                 {
//                     label: 'General',
//                     children: [
//                         { label: 'Nautical', layer: usa_naut },
//                         { label: 'Satellite', layer: usa_sat },
//                         { label: 'Topographical', layer: usa_topo },
//                     ]
//                 },
//                 {
//                     label: 'States',
//                     children: [
//                         { label: 'CA', layer: usa_ca },
//                         { label: 'NY', layer: usa_ny },
//                         /* ... */
//                     ]
//                 }
//             ]
//         },
//     ]
// };

// var overlaysTree = {
//     label: 'Points of Interest',
//     selectAllCheckbox: 'Un/select all',
//     children: [
//         {
//             label: 'Europe',
//             selectAllCheckbox: true,
//             children: [
//                 {
//                     label: 'France',
//                     selectAllCheckbox: true,
//                     children: [
//                         { label: 'Tour Eiffel', layer: L.marker([27000, 27000]) },
//                         { label: 'Notre Dame', layer: L.marker([47000, 7000]) },
//                         { label: 'Louvre', layer: L.marker([7000, 1000]) },
//                     ]
//                 }, {
//                     label: 'Germany',
//                     selectAllCheckbox: true,
//                     children: [
//                         { label: 'Branderburger Tor', layer: L.marker([0, 0])},
//                         { label: 'Kölner Dom', layer: L.marker([50000, 60000])},
//                     ]
//                 }, {label: 'Spain',
//                     selectAllCheckbox: 'De/seleccionar todo',
//                     children: [
//                         { label: 'Palacio Real', layer: L.marker([90000, 27000])},
//                         { label: 'La Alhambra', layer: L.marker([27000, 100000])},
//                     ]
//                 }
//             ]
//         }, {
//             label: 'Asia',
//             selectAllCheckbox: true,
//             children: [
//                 {
//                     label: 'Jordan',
//                     selectAllCheckbox: true,
//                     children: [
//                         { label: 'Petra', layer: L.marker([95000, 100000]) },
//                         { label: 'Wadi Rum', layer: L.marker([45000, 11000]) }
//                     ]
//                 }, {
//                 /* ... */
//                 }
//             ]
//         }
//     ]
// }

// L.control.layers.tree({}, overlaysTree, {position:'topleft'}).addTo(map);