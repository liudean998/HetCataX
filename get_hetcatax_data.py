import requests


base_url = "https://hetcatax.com"


def get_token():
    user_name = '***'
    password = '***'
    login_url = f'{base_url}/api/signin/'
    headers = {'Email': user_name, 'Password': password}
    response = requests.post(login_url, headers=headers)
    token = response.headers.get('Authorization')
    return token


def get_data():
    token = get_token()
    headers = {'Authorization': token}

    payload = {
        'base': {
            'model': 1,  # 查询模式 0： 至少包含， 1：只包含，2：使用formula查询
            'element': ['Cu', 'Zn'], # 设置目标元素
            'formula': {}  # 精确配比，比如{"Cu": 24, "Y": 2}，使用这个信息查询需要model设置为2
        },
        "abstract": {
            "num_of_atoms": [20, 100],  # 原子数上下限
            "source": "" # 数据来源
        },
        "model_type": "catalyst", # 数据模块
        "size": 20
    }

    response = requests.post(f'{base_url}/api/self_g/r_relax/', json=payload, headers=headers)
    data = response.json()['data']['infos']
    for i in data:
        print(i)
    print(data)
    print(len(data))
    return data


get_data()
