export const PLOT_CONFIG = {
    margin: {
        top: 60,
        right: 80,
        bottom: 50,
        left: 100
    },
    point: {
        radius: 5,
        // color: 'blue',
        color: '#69b3a2',
        hoverColor: '#4a90e2'
    },
    diagonalLine: {
        color: 'red',
        width: 2,
        dashArray: '5,5'
    },
    animation: {
        duration: 1000,
        tooltip: {
            fadeIn: 200,
            fadeOut: 500
        }
    },
    guide: {
        backgroundColor: 'rgba(255, 223, 186, 0.7)',
        padding: 5,
        lineHeight: 15,
        fontSize: '12px'
    },
    axis: {
        fontSize: '12px',
        tickSize: 5,
        tickPadding: 5
    }
};

export function calculateDimensions() {
    return {
        width: window.innerWidth - PLOT_CONFIG.margin.left - PLOT_CONFIG.margin.right,
        height: window.innerHeight * 0.34 - PLOT_CONFIG.margin.top - PLOT_CONFIG.margin.bottom
    };
}
