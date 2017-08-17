#!/bin/env python3

## Tommy Carstensen, Louise Barr, July 2017
## Genome Research Ltd.

## Import the requests module to communicate with FitBit.
import requests
## Import the sys module to parse command line options.
## Module argparse also an option.
import sys
## Import the module os to check for file existence.
import os
## Import the module pandas to read csv files and create data frames.
import pandas as pd
## Import datetime to convert date string to datetime date object.
import datetime
## Import matplotlib to make plots.
import matplotlib.pyplot as plt


def main():

    ## Parse command line options.
    consumer_key = sys.argv[1]
    assert len(consumer_key) == 6
    date1 = sys.argv[2]
    assert len(date1) == 10
    access_token = sys.argv[3]
    sangerID = sys.argv[4]

    date2 = datetime.datetime.strftime(datetime.datetime.strptime(
        date1, '%Y-%m-%d') + datetime.timedelta(days=1), '%Y-%m-%d')

    ## https://dev.fitbit.com/docs/activity/#get-activity-intraday-time-series
    parse_fitbit('steps', date1, date2, consumer_key, access_token, '1min')

    ## https://dev.fitbit.com/docs/heart-rate/#get-heart-rate-intraday-time-series
    parse_fitbit('heart', date1, date2, consumer_key, access_token, '1sec')

    ## Read Actigraph data into a Pandas dataframe.
    with open('{}.csv'.format(sangerID)) as f:
        for line in f:
            if line.startswith('Start Time'):
                break
        start_time = line.rstrip().split()[2]
    df = pd.read_csv(
        '{}.csv'.format(sangerID),
        header=10,
        usecols=('Steps', 'HR',),
        )
##    df.index = pd.date_range(date+' '+start_time, periods=len(df), freq='s')
    df.index = pd.to_datetime(df.index, origin=date1+' '+start_time, unit='s')

    df1 = parse_csv(
        '{}.{}.heart.csv'.format(consumer_key, date1), ('HR_Fitbit',), date1,)
    df2 = parse_csv(
        '{}.{}.heart.csv'.format(consumer_key, date2), ('HR_Fitbit',), date2,)
    df_1sec = df.join(pd.concat([df1, df2]), how='inner')
    df_1sec = df_1sec[df_1sec.HR > 0]  # Tell Louise about this!!! Actigraph many zero!!!
    del df_1sec['Steps']

    df1 = parse_csv(
        '{}.{}.steps.csv'.format(consumer_key, date1), ('Steps_Fitbit',), date1,)
    df2 = parse_csv(
        '{}.{}.steps.csv'.format(consumer_key, date2), ('Steps_Fitbit',), date2,)
    df_1min = df.resample('T').sum().join(pd.concat([df1, df2]), how='inner')
    df_1min = df_1min[(df_1min.Steps > 0) & (df_1min.Steps_Fitbit > 0)]
    del df_1min['HR']

    print('raw', df.tail())
    print('1sec', df_1sec.tail())
    print('1min', df_1min.tail())

    plots(df_1min, df_1sec, consumer_key, sangerID, date1)

    return


def plots(df_1min, df_1sec, consumer_key, sangerID, date1):

    for df, activity in ((df_1min, 'Steps'), (df_1sec, 'HR')):
        ## Calculate Pearson correlation coefficient.
        corr = pd.DataFrame.corr(df, method='pearson')
        ## Create title for all plots.
        title = '{}, {}, {}, {}'.format(activity, consumer_key, sangerID, date1)
        ## Plot time series.
        plt.clf()
        fig, ax = plt.subplots()
        plot = df[activity].plot()
        plot = df[activity+'_Fitbit'].plot()
    ##    fig = plot.get_figure()
        plt.title(title)
        plt.legend()
        ax.text(
            0.01, 0.99,
            'n = {}\nr = {:.2f}'.format(
                len(df), corr[activity][activity+'_Fitbit']),
            verticalalignment='top', horizontalalignment='left',
            transform=ax.transAxes,
            color='green', fontsize='small',
            )
        plt.savefig('{}.timeseries.{}.png'.format(consumer_key, activity))
        ## Calculate mean and diff for BA plot.
        df['mean'] = df.mean(axis=1)
        df['diff'] = df[activity] - df[activity+'_Fitbit']
        ## Plot Bland-Altman / Tukey Mean Difference plots.
        plt.clf()
        ax = df.plot('mean', 'diff', kind='scatter')
        ax.set_title(title)
        ax.set_xlabel('mean')
        ax.set_ylabel('diff')
        ax.text(
            0.99, 0.01,
            'n = {}\nr = {:.2f}'.format(
                len(df), corr[activity][activity+'_Fitbit']),
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes,
            color='green', fontsize='small',
            )
        plt.savefig('{}.BA.{}.png'.format(consumer_key, activity))
        ## Plot scatter plot.
        plt.clf()
        ax = df.plot(activity, activity+'_Fitbit', kind='scatter')
        ax.set_title(title)
        ax.set_xlabel(activity)
        ax.set_ylabel(activity+'_Fitbit')
        ax.text(
            0.99, 0.01,
            'n = {}\nr = {:.2f}'.format(
                len(df), corr[activity][activity+'_Fitbit']),
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes,
            color='green', fontsize='small',
            )
        plt.savefig('{}.scatter.{}.png'.format(consumer_key, activity))

    return


def parse_fitbit(activity, date1, date2, consumer_key, access_token, detail_level):

    resource_path = 'activities/{activity}'.format(activity=activity)

    for date in (date1, date2):
        url = 'https://api.fitbit.com/1/user/{consumer_key}/{resource_path}/date/{date}/1d/{detail_level}.json'.format(
            date=date1, consumer_key=consumer_key, resource_path=resource_path,
            end_date=date, detail_level=detail_level,
            )
        if os.path.isfile('{}.{}.{}.csv'.format(consumer_key, date, activity)):
            continue
        response = requests.get(
            url=url, headers={'Authorization': 'Bearer ' + access_token})
        if not response.ok:
            print_response_and_exit(response)
        print(date, activity, response.json()['activities-'+activity][0]['value'])
        dataset = response.json()['activities-{}-intraday'.format(activity)]['dataset']
        with open('{}.{}.{}.csv'.format(consumer_key, date, activity), 'w') as f:
            for d in dataset:
                print(d['time'], d['value'], sep=',', file=f)

    return


def parse_csv(path_csv, names, date):

    df = pd.read_csv(
        path_csv,
        header=None,
        index_col=0,
        names=names,
        )
    df.index = pd.to_datetime(date+' '+df.index, format='%Y-%m-%d %H:%M:%S')

    return df


def print_response_and_exit(response):

    print(response)
    print('ok', response.ok)
    print('status_code', response.status_code)
    print('reason', response.reason)
    exit()

    return

if __name__ == '__main__':
    main()

## 1) Register an app
## https://dev.fitbit.com/apps/new
## 2) Manage my apps (i.e. get Client ID and Client Secret)
## https://dev.fitbit.com/apps
## 3a) OAuth 2.0 tutorial page
##  1 Authorize (redirect user to authorization page)
##  1A Get Code (after redirection from authroization page back to redirect URI)
##  Parse response
##  Make Request
##  Refresh Token (when needed)
## 3b) OR while ensuring Callback URL of app is http://127.0.0.1:8080/
##  feed app client_id and client_secret
##  to https://github.com/orcasgit/python-fitbit/blob/master/gather_keys_oauth2.py
## 3c) OR https://github.com/niwasawa/fitbit-api-python-client/blob/master/get_access_token.py

## https://python-fitbit.readthedocs.io/en/latest/
## https://github.com/orcasgit/python-fitbit

## --------------------------------------------------

## https://dev.fitbit.com/docs/basics/#personal
## "Personal" applications authenticate using either the Authorization Code Grant Flow or the Implicit Grant Flow, and also have access to the intraday time series data (such as Heart Rate or Activity) of the owner of the app only.

## https://dev.fitbit.com/docs/oauth2/#obtaining-consent
## Applications that do not have a web service should use the Implicit Grant flow.

## https://dev.fitbit.com/docs/oauth2/#implicit-grant-flow
## Unlike the Authorization Code Grant Flow, the refresh tokens are not issued with the Implicit Grant flow. Refreshing a token requires use of the client secret, which cannot safely be stored in distributed application code. When the access token expires, users will need to re-authorize your app.
## Access tokens obtained via the Implicit Grant Flow should only be stored on the device used to obtain the authorization. If your application has a web server component, your application should use the Authorization Code Grant Flow.

## https://dev.fitbit.com/docs/oauth2/#authorization-page
## Examples
## Authorization Code Flow:
## https://www.fitbit.com/oauth2/authorize?response_type=code&client_id=22942C&redirect_uri=http%3A%2F%2Fexample.com%2Ffitbit_auth&scope=activity%20nutrition%20heartrate%20location%20nutrition%20profile%20settings%20sleep%20social%20weight
## Implicit Grant Flow:
## https://www.fitbit.com/oauth2/authorize?response_type=token&client_id=22942C&redirect_uri=http%3A%2F%2Fexample.com%2Ffitbit_auth&scope=activity%20nutrition%20heartrate%20location%20nutrition%20profile%20settings%20sleep%20social%20weight&expires_in=604800

## https://dev.fitbit.com/docs/oauth2/#refreshing-tokens
## When an access token expires, an HTTP 401 error will be returned:
## Your application will need to refresh the access token. The Fitbit API follows the OAuth 2.0 specification for refreshing access tokens. A refresh token does not expire until it is used. A refresh token can only be used once, as a new refresh token is returned with the new access token.

## https://dev.fitbit.com/docs/activity/#get-activity-intraday-time-series
## Access to the Intraday Time Series for personal use (accessing your own data) is available through the "Personal" App Type.
## Access to the Intraday Time Series for all other uses is currently granted on a case-by-case basis. Applications must demonstrate necessity to create a great user experience. Fitbit is very supportive of non-profit research and personal projects.
