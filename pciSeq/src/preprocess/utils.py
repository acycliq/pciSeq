
def get_img_shape(coo):
    img_shape = set([d.shape for d in coo])
    assert len(img_shape) == 1, 'pages do not have the same shape'
    img_shape = img_shape.pop()
    w = img_shape[1]
    h = img_shape[0]
    return [h, w]
