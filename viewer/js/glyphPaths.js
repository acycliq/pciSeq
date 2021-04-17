
function ctxPath(glyphName, ctx, p, r) {

    if (glyphName === 'star6') {
        ctx.beginPath();
        ctx.moveTo(p.x + r, p.y);
        ctx.lineTo(p.x + 0.43 * r, p.y + 0.25 * r);
        ctx.lineTo(p.x + 0.50 * r, p.y + 0.87 * r);
        ctx.lineTo(p.x, p.y + 0.50 * r);
        ctx.lineTo(p.x - 0.50 * r, p.y + 0.87 * r);
        ctx.lineTo(p.x - 0.43 * r, p.y + 0.25 * r);
        ctx.lineTo(p.x - r, p.y);
        ctx.lineTo(p.x - 0.43 * r, p.y - 0.25 * r);
        ctx.lineTo(p.x - 0.50 * r, p.y - 0.87 * r);
        ctx.lineTo(p.x, p.y - 0.50 * r);
        ctx.lineTo(p.x + 0.50 * r, p.y - 0.87 * r);
        ctx.lineTo(p.x + 0.43 * r, p.y - 0.25 * r);
        ctx.closePath();
    }


    else if (glyphName === 'star5') {
        ctx.beginPath();
        ctx.moveTo(p.x, p.y + 0.658351875 * r);
        ctx.lineTo(p.x + 0.618027443 * r, p.y + 1 * r);
        ctx.lineTo(p.x + 0.5 * r, p.y + 0.27637816 * r);
        ctx.lineTo(p.x + 1 * r, p.y - 0.236080209 * r);
        ctx.lineTo(p.x + 0.309026865 * r, p.y - 0.341634306 * r);
        ctx.lineTo(p.x + 0 * r, p.y - 1 * r);
        ctx.lineTo(p.x - 0.309026865 * r, p.y - 0.341634306 * r);
        ctx.lineTo(p.x - 1 * r, p.y - 0.236080209 * r);
        ctx.lineTo(p.x - 0.5 * r, p.y + 0.27637816 * r);
        ctx.lineTo(p.x - 0.618027443 * r, p.y + 1 * r);
        ctx.lineTo(p.x, p.y + 0.658351875 * r);
        ctx.closePath();
    }


    else if (glyphName === 'diamond') {
        ctx.beginPath();
        ctx.moveTo(p.x - r, p.y);
        ctx.lineTo(p.x, p.y - r);
        ctx.lineTo(p.x + r, p.y);
        ctx.lineTo(p.x, p.y + r);
        ctx.lineTo(p.x - r, p.y);
        ctx.closePath();
    }


    else if (glyphName === 'square') {
        ctx.beginPath();
        ctx.moveTo(p.x - r, p.y - r);
        ctx.lineTo(p.x + r, p.y - r);
        ctx.lineTo(p.x + r, p.y + r);
        ctx.lineTo(p.x - r, p.y + r);
        ctx.lineTo(p.x - r, p.y - r);
        ctx.closePath();
    }


    else if (glyphName === 'triangleUp') {
        ctx.beginPath();
        ctx.moveTo(p.x - r, p.y + r);
        ctx.lineTo(p.x, p.y - r);
        ctx.lineTo(p.x + r, p.y + r);
        ctx.closePath();
    }


    else if (glyphName === 'triangleDown') {
        ctx.beginPath();
        ctx.moveTo(p.x - r, p.y - r);
        ctx.lineTo(p.x, p.y + r);
        ctx.lineTo(p.x + r, p.y - r);
        ctx.closePath();
    }


    else if (glyphName === 'triangleRight') {
        ctx.beginPath();
        ctx.moveTo(p.x - r, p.y - r);
        ctx.lineTo(p.x + r, p.y);
        ctx.lineTo(p.x - r, p.y + r);
        ctx.closePath();
    }


    else if (glyphName === 'triangleLeft') {
        ctx.beginPath();
        ctx.moveTo(p.x + r, p.y - r);
        ctx.lineTo(p.x - r, p.y);
        ctx.lineTo(p.x + r, p.y + r);
        ctx.closePath();
    }


    else if (glyphName === 'cross') {
        ctx.beginPath();
        ctx.moveTo(p.x + r, p.y + r);
        ctx.lineTo(p.x - r, p.y - r);
        ctx.moveTo(p.x - r, p.y + r);
        ctx.lineTo(p.x + r, p.y - r);
        ctx.closePath();
    }


    else if (glyphName === 'plus') {
        ctx.beginPath();
        ctx.moveTo(p.x, p.y + r);
        ctx.lineTo(p.x, p.y - r);
        ctx.moveTo(p.x - r, p.y);
        ctx.lineTo(p.x + r, p.y);
        ctx.closePath();
    }


    else if (glyphName === 'asterisk') {
        ctx.beginPath();
        ctx.moveTo(p.x, p.y + r);
        ctx.lineTo(p.x, p.y - r);
        ctx.moveTo(p.x - r, p.y);
        ctx.lineTo(p.x + r, p.y);
        ctx.moveTo(p.x + 0.5 * r, p.y + 0.5 * r);
        ctx.lineTo(p.x - 0.5 * r, p.y - 0.5 * r);
        ctx.moveTo(p.x - 0.5 * r, p.y + 0.5 * r);
        ctx.lineTo(p.x + 0.5 * r, p.y - 0.5 * r);
        ctx.arc(p.x, p.y, 2, 0, Math.PI * 2, true);
        ctx.closePath();
    }


    else if (glyphName === 'circle') {
        ctx.beginPath();
        ctx.arc(p.x, p.y, r, 0, Math.PI * 2, true);
        ctx.closePath();
    }


    else if (glyphName === 'point') {
        ctx.beginPath();
        ctx.arc(p.x, p.y, 3, 0, Math.PI * 2, true);
        ctx.arc(p.x, p.y, 2, 0, Math.PI * 2, true);
        ctx.closePath();
    }

    else {
        console.log('glyph: "' + glyphName + '" not implemented.')
    }

    return ctx
}




