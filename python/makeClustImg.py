import numpy as np
import matplotlib.pyplot as plt


def make_image(pix_array):
    image_array_layer1_single_phi = np.zeros((1920, 4992))
    image_array_layer4_single_phi = np.zeros((10240, 4992))
    image_array_layer1_double_phi = np.zeros((3840, 4992))
    image_array_layer4_double_phi = np.zeros((20480, 4992))

    image_array_layer1_single_phi_double_z = np.zeros((1920, 9984))
    image_array_layer4_single_phi_double_z = np.zeros((10240, 9984))
    image_array_layer1_double_phi_double_z = np.zeros((3840, 9984))
    image_array_layer4_double_phi_double_z = np.zeros((20480, 9984))

    for pix in pix_array[0]:
        phi = pix[0]  # singles 0.00327, 0.000614; doubles 0.00164, 0.000307
        z = pix[1]  # singles 0.0168; doubles 0.00841
        r = pix[2]

        phi = phi + np.pi
        z = z + 28

        phi_single_layer1 = int(phi / (2*np.pi/1920.))
        phi_double_layer1 = int(phi / (2*np.pi/3840.))
        phi_single_layer4 = int(phi / (2*np.pi/10240.))
        phi_double_layer4 = int(phi / (2*np.pi/20480.))

        z_single = int(z / (56./4992.))
        z_double = int(z / (56./9984.))

        image_array_layer1_single_phi[phi_single_layer1][z_single] = r
        image_array_layer1_double_phi[phi_double_layer1][z_single] = r
        image_array_layer4_single_phi[phi_single_layer4][z_single] = r
        image_array_layer4_double_phi[phi_double_layer4][z_single] = r

        image_array_layer1_single_phi_double_z[phi_single_layer1][z_double] = r
        image_array_layer1_double_phi_double_z[phi_double_layer1][z_double] = r
        image_array_layer4_single_phi_double_z[phi_single_layer4][z_double] = r
        image_array_layer4_double_phi_double_z[phi_double_layer4][z_double] = r

    plt.imshow(image_array_layer1_single_phi, cmap=plt.cm.bone)
    plt.savefig('image_layer1_single_phi.jpg', bbox_inches='tight')
    plt.close()

    plt.imshow(image_array_layer1_double_phi, cmap=plt.cm.bone)
    plt.savefig('image_layer1_double_phi.jpg', bbox_inches='tight')
    plt.close()

    plt.imshow(image_array_layer4_single_phi, cmap=plt.cm.bone)
    plt.savefig('image_layer4_single_phi.jpg', bbox_inches='tight')
    plt.close()

    plt.imshow(image_array_layer4_double_phi, cmap=plt.cm.bone)
    plt.savefig('image_layer4_double_phi.jpg', bbox_inches='tight')
    plt.close()

    plt.imshow(image_array_layer1_single_phi_double_z, cmap=plt.cm.bone)
    plt.savefig('image_layer1_single_phi_double_z.jpg', bbox_inches='tight')
    plt.close()

    plt.imshow(image_array_layer1_double_phi_double_z, cmap=plt.cm.bone)
    plt.savefig('image_layer1_double_phi_double_z.jpg', bbox_inches='tight')
    plt.close()

    plt.imshow(image_array_layer4_single_phi_double_z, cmap=plt.cm.bone)
    plt.savefig('image_layer4_single_phi_double_z.jpg', bbox_inches='tight')
    plt.close()

    plt.imshow(image_array_layer4_double_phi_double_z, cmap=plt.cm.bone)
    plt.savefig('image_layer4_double_phi_double_z.jpg', bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    in_array = np.load('phi_z_r_array.npy')
    make_image(in_array)
