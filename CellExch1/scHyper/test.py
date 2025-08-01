import pickle
import pkg_resources
file_path = pkg_resources.resource_filename(
    __name__, f"predictions.pkl")
with open(file_path, 'rb') as file:
    predictions = pickle.load(file)
predictions = predictions.cpu()
#df_nn = pd.DataFrame(predictions.numpy(), columns=['Column1'])
print(predictions)
print(predictions.shape)