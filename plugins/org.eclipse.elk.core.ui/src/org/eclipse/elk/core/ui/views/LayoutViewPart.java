/*******************************************************************************
 * Copyright (c) 2009, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.ui.views;

import java.util.LinkedList;
import java.util.List;

import org.eclipse.core.runtime.IStatus;
import org.eclipse.core.runtime.Status;
import org.eclipse.elk.core.service.LayoutConnectorsService;
import org.eclipse.elk.core.ui.ElkUiPlugin;
import org.eclipse.elk.core.util.Maybe;
import org.eclipse.jface.action.ActionContributionItem;
import org.eclipse.jface.action.IToolBarManager;
import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.jface.viewers.ISelection;
import org.eclipse.swt.layout.FormAttachment;
import org.eclipse.swt.layout.FormData;
import org.eclipse.swt.layout.FormLayout;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.swt.widgets.Control;
import org.eclipse.swt.widgets.Display;
import org.eclipse.swt.widgets.Tree;
import org.eclipse.swt.widgets.TreeItem;
import org.eclipse.ui.IActionBars;
import org.eclipse.ui.IEditorPart;
import org.eclipse.ui.IPartListener;
import org.eclipse.ui.ISelectionListener;
import org.eclipse.ui.IWorkbenchPage;
import org.eclipse.ui.IWorkbenchPart;
import org.eclipse.ui.IWorkbenchWindow;
import org.eclipse.ui.PlatformUI;
import org.eclipse.ui.forms.widgets.Form;
import org.eclipse.ui.forms.widgets.FormToolkit;
import org.eclipse.ui.part.ViewPart;
import org.eclipse.ui.statushandlers.StatusManager;
import org.eclipse.ui.views.properties.IPropertySheetEntry;
import org.eclipse.ui.views.properties.PropertySheetPage;

import com.google.inject.ConfigurationException;
import com.google.inject.Injector;

/**
 * A view that displays layout options for selected objects.
 *
 * @author msp
 */
public class LayoutViewPart extends ViewPart {
    
    /**
     * Finds the active layout view, if it exists.
     * 
     * @return the active layout view, or {@code null} if there is none
     */
    public static LayoutViewPart findView() {
        if (Display.getCurrent() == null) {
            final Maybe<LayoutViewPart> part = new Maybe<LayoutViewPart>();
            Display.getDefault().syncExec(new Runnable() {
                public void run() {
                    part.set(findViewUI());
                }
            });
            return part.get();
        } else {
            return findViewUI();
        }
    }
    
    /**
     * Finds the active layout view, if it exists. This method works only if called
     * from the user interface thread.
     * 
     * @return the active layout view, or {@code null} if there is none
     */
    private static LayoutViewPart findViewUI() {
        IWorkbenchWindow activeWindow = PlatformUI.getWorkbench().getActiveWorkbenchWindow();
        if (activeWindow != null) {
            IWorkbenchPage activePage = activeWindow.getActivePage();
            if (activePage != null) {
                return (LayoutViewPart) activePage.findView(VIEW_ID);
            }
        }
        return null;
    }

    /** The view identifier. */
    public static final String VIEW_ID = "org.eclipse.elk.core.views.layout";
    /** Preference identifier for enabling categories. */
    public static final String PREF_CATEGORIES = "view.categories";
    /** Preference identifier for enabling advanced options. */
    public static final String PREF_ADVANCED = "view.advanced";
    
    /** The form toolkit used to create forms. */
    private FormToolkit toolkit;
    /** The page that is displayed in this view part. */
    private PropertySheetPage page;
    /** The last created property source provider. */
    private LayoutPropertySourceProvider currentPropSourceProvider;
    
    /** The selection listener. */
    private final ISelectionListener selectionListener = (final IWorkbenchPart part, final ISelection selection) -> {
        Injector injector = LayoutConnectorsService.getInstance().getInjector(part, null);
        if (injector != null) {
            try {
                currentPropSourceProvider = injector.getInstance(LayoutPropertySourceProvider.class);
                currentPropSourceProvider.setWorkbenchPart(part);
                page.setPropertySourceProvider(currentPropSourceProvider);
                page.selectionChanged(part, selection);
            } catch (ConfigurationException exception) {
                IStatus status = new Status(IStatus.ERROR, ElkUiPlugin.PLUGIN_ID,
                        "The Guice configuration for " + part.getTitle() + " is inconsistent.",
                        exception);
                StatusManager.getManager().handle(status, StatusManager.LOG);
            }
        } else if (part instanceof IEditorPart) {
            currentPropSourceProvider = null;
            page.setPropertySourceProvider(null);
        }
    };
    
    /** The part listener for reacting to closed workbench parts. */
    private final IPartListener partListener = new IPartListener() {
        public void partClosed(final IWorkbenchPart part) {
            if (currentPropSourceProvider != null && part == currentPropSourceProvider.getWorkbenchPart()) {
                // Reset the property source provider in order to ensure proper garbage collection
                currentPropSourceProvider = null;
                page.setPropertySourceProvider(null);
            }
        }
        public void partOpened(final IWorkbenchPart part) { }
        public void partDeactivated(final IWorkbenchPart part) { }
        public void partBroughtToTop(final IWorkbenchPart part) { }
        public void partActivated(final IWorkbenchPart part) { }
    };
    
    /**
     * Return the currently active property source provider, or {@code null} if none is active.
     */
    public LayoutPropertySourceProvider getPropertySourceProvider() {
        return currentPropSourceProvider;
    }
    
    /** margin width for the form layout. */
    private static final int MARGIN_WIDTH = 4;
    /** position for left attachment. */
    private static final int FORM_LEFT = 0;
    /** position for right attachment. */
    private static final int FORM_RIGHT = 100;
    /** position for top attachment. */
    private static final int FORM_TOP = 0;
    /** position for bottom attachment. */
    private static final int FORM_BOTTOM = 100;
    
    @Override
    public void createPartControl(final Composite parent) {
        // GUI code needs magic numbers
        // CHECKSTYLEOFF MagicNumber
        
        toolkit = new FormToolkit(parent.getDisplay());
        Form form = toolkit.createForm(parent);
        
        // Content
        Composite content = form.getBody();
        FormLayout contentLayout = new FormLayout();
        contentLayout.marginWidth = MARGIN_WIDTH;
        content.setLayout(contentLayout);
        
        // Property Sheet Page
        page = new PropertySheetPage();
        page.setRootEntry(new ValidatingPropertySheetEntry());
        page.createControl(content);
        FormData formData = new FormData();
        formData.left = new FormAttachment(FORM_LEFT, 0);
        formData.right = new FormAttachment(FORM_RIGHT, 0);
        formData.top = new FormAttachment(FORM_TOP, 5);
        formData.bottom = new FormAttachment(FORM_BOTTOM, 0);
        page.getControl().setLayoutData(formData);
        IPreferenceStore preferenceStore = ElkUiPlugin.getInstance().getPreferenceStore();
        
        // add actions to the toolbar, view menu, and context menu
        IActionBars actionBars = getViewSite().getActionBars();
        page.setActionBars(actionBars);
        IToolBarManager toolBarManager = actionBars.getToolBarManager();
        
        // CHECKSTYLEON MagicNumber
        // set the stored value of the categories button
        ActionContributionItem categoriesItem =
                (ActionContributionItem) toolBarManager.find("categories");
        if (categoriesItem != null) {
            categoriesItem.getAction().setChecked(preferenceStore.getBoolean(PREF_CATEGORIES));
            categoriesItem.getAction().run();
        }
        
        // set the stored value of the advanced button
        ActionContributionItem advancedItem =
                (ActionContributionItem) toolBarManager.find("filter");
        if (advancedItem != null) {
            advancedItem.getAction().setChecked(preferenceStore.getBoolean(PREF_ADVANCED));
            advancedItem.getAction().run();
        }
        
        // get the current selection and trigger an update
        IWorkbenchWindow workbenchWindow = getSite().getWorkbenchWindow();
        IWorkbenchPage activePage = workbenchWindow.getActivePage();
        if (activePage != null) {
            final IWorkbenchPart activePart = activePage.getActivePart();
            final ISelection selection = workbenchWindow.getSelectionService().getSelection();
            if (activePart != null && selection != null) {
                workbenchWindow.getWorkbench().getDisplay().asyncExec(new Runnable() {
                    public void run() {
                        selectionListener.selectionChanged(activePart, selection);
                    }
                });
            }
        }
        workbenchWindow.getSelectionService().addSelectionListener(selectionListener);
        workbenchWindow.getPartService().addPartListener(partListener);
    }

    @Override
    public void setFocus() {
        page.setFocus();
    }
    
    /**
     * Returns the control that is handled by this view part.
     * 
     * @return the control
     */
    public Control getControl() {
        return page.getControl();
    }
    
    @Override
    public void dispose() {
        // Store the current status of the categories button
        ActionContributionItem categoriesItem = (ActionContributionItem) getViewSite()
                .getActionBars().getToolBarManager().find("categories");
        if (categoriesItem != null) {
            ElkUiPlugin.getInstance().getPreferenceStore().setValue(PREF_CATEGORIES,
                    categoriesItem.getAction().isChecked());
        }
        // Store the current status of the advanced button
        ActionContributionItem advancedItem = (ActionContributionItem) getViewSite()
                .getActionBars().getToolBarManager().find("filter");
        if (advancedItem != null) {
            ElkUiPlugin.getInstance().getPreferenceStore().setValue(PREF_ADVANCED,
                    advancedItem.getAction().isChecked());
        }
        // Dispose the view part
        getSite().getWorkbenchWindow().getSelectionService().removeSelectionListener(selectionListener);
        getSite().getWorkbenchWindow().getPartService().removePartListener(partListener);
        currentPropSourceProvider = null;
        toolkit.dispose();
        super.dispose();
    }
    
    /**
     * Refreshes the layout view asynchronously.
     */
    public void refresh() {
        Display.getDefault().asyncExec(new Runnable() {
            public void run() {
                page.refresh();
            }
        });
    }
    
    /**
     * Returns the current selection of the layout view.
     * 
     * @return the selected property sheet entries
     */
    public List<IPropertySheetEntry> getSelection() {
        List<IPropertySheetEntry> entries = new LinkedList<IPropertySheetEntry>();
        TreeItem[] treeItems = ((Tree) page.getControl()).getSelection();
        for (TreeItem item : treeItems) {
            Object data = item.getData();
            if (data instanceof IPropertySheetEntry) {
                entries.add((IPropertySheetEntry) data);
            } else {
                // A category was selected, apply options for all children
                for (TreeItem childItem : item.getItems()) {
                    data = childItem.getData();
                    if (data instanceof IPropertySheetEntry) {
                        entries.add((IPropertySheetEntry) data);
                    }
                }
            }
        }
        return entries;
    }

}
